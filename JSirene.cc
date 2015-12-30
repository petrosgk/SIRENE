#include <string>
#include <iostream>
#include <iomanip>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TProfile.h"
#include "TRandom3.h"

#include "antcc/Header.hh"
#include "antcc/Event.hh" 

#include "JROOT/JRootTree.hh"

#include "JLang/JSharedPointer.hh"
#include "Jeep/JTimer.hh"
#include "Jeep/JParser.hh"
#include "Jeep/JMessage.hh"
#include "JPhysics/JPDFToolkit.hh"
#include "JPhysics/JCDFTable.hh"
#include "JPhysics/JRadiationSource.hh"
#include "JTools/JFunction1D_t.hh"
#include "JTools/JFunctionalMap_t.hh"

#include "JSirene/JSirene.hh"
#include "JSirene/JEvent.hh"
#include "JSirene/JPythia.hh"
#include "JSirene/JSeaWater.hh"
#include "JSirene/JCDFTable1D.hh"
#include "JSirene/JCDFTable2D.hh"
#include "JDetector/JDetector.hh"
#include "JDetector/JDetectorSubset.hh"
#include "JDetector/JDetectorToolkit.hh"
#include "JGeometry3D/JCylinder3D.hh"

#include "JSupport/JMultipleFileScanner.hh"
#include "JSupport/JFileRecorder.hh"
#include "JSupport/JMonteCarloToolkit.hh"
#include "JSupport/JSupport.hh"

#include "omp.h"

/**   --------------------------------------------------------------------------
 *   |                                                                          |
 *   |  Program to simulate detector response to muons and showers.             |
 *   |                                                                          |
 *   |   - CDF file descriptor should contain '%' character;                    |
 *   |     The file names are obtained by replacing % with 1, 2, 13 and 14;     |
 *   |                                                                          |
 *   |   - Detector geometry is read from .det file;                            |
 *   |                                                                          |
 *   |   - Events are read from / written to ROOT formatted .evt file;          |
 *    --------------------------------------------------------------------------
 */
int main(int argc, char **argv) {
	using namespace std;
	using namespace JSUPPORT;

	string fileDescriptor;
	JMultipleFileScanner<JMonteCarloTypes_t> inputFile;
	JFileRecorder<JLANG::JAppend<JMonteCarloTypes_t, TRandom>::typelist> outputFile;
	JLimit_t& numberOfEvents = inputFile.getLimit();
	string detectorFile;
	double Tmax;
	bool geasim;
	bool writeEMShowers;
	UInt_t seed;
	int debug;

	// --- enable OpenMP ---
	bool OMP = true;
	if (!OMP) {
		omp_set_num_threads(1);
	} else {
		omp_set_num_threads(4);
	}
	// --- enable OpenMP ---

	try {

		JParser<> zap;

		zap['F'] = make_field(fileDescriptor);
		zap['f'] = make_field(inputFile);
		zap['o'] = make_field(outputFile)= "sirene.root";
		zap['n'] = make_field(numberOfEvents)= JLimit::max();
		zap['a'] = make_field(detectorFile);
		zap['T'] = make_field(Tmax)= 0.0;
		zap['G'] = make_field(geasim);
		zap['s'] = make_field(writeEMShowers);
		zap['S'] = make_field(seed)= 0;
		zap['d'] = make_field(debug)= 1;

		if (zap.read(argc, argv) != 0)
			return 1;
	} catch (exception &error) {
		FATAL(error.what() << endl);
	}

	gRandom->SetSeed(seed);

	using namespace JTOOLS;
	using namespace JPHYSICS;
	using namespace JSIRENE;
	using namespace JEEP;
	using namespace JLANG;
	using namespace JGEOMETRY3D;

	typedef JSplineFunction1D_t JFunction1D_t;
	typedef JMapList<JPolint1FunctionalMap,
			JMapList<JPolint1FunctionalGridMap,
					JMapList<JPolint1FunctionalGridMap> > > J3DMap_t;
	typedef JMapList<JPolint1FunctionalMap,
			JMapList<JPolint1FunctionalMap,
					JMapList<JPolint1FunctionalGridMap,
							JMapList<JPolint1FunctionalGridMap> > > > J4DMap_t;

	typedef JCDFTable<JFunction1D_t, J3DMap_t> JCDF4D_t;
	typedef JCDFTable<JFunction1D_t, J4DMap_t> JCDF5D_t;

	JCDF4D_t f1;       // direct    light from muon
	JCDF4D_t f2;       // scattered light from muon
	JCDF5D_t f13;      // direct    light from EM shower
	JCDF5D_t f14;      // scattered light from EM shower

	double Ecut = 0.1; // minimal energy for generation of light from shower [GeV]
	double Emin = 1.0; // minimal energy of muon for shower generation [GeV]
	double Zbed = 0.0;                             // level of seabed [m]

	for (string buffer[] = { "1", "2", "13", "14", "" }, *k = buffer; *k != "";
			++k) {

		string file_name = fileDescriptor;

		string::size_type ipos = file_name.find('%');

		if (ipos == string::npos) {
			FATAL("error file descriptor " << fileDescriptor << endl);
		}

		file_name.replace(ipos, 1, *k);

		NOTICE("loading input from file " << file_name << "... " << flush);

		try {
			if (*k == "1")
				f1.load(file_name.c_str());
			if (*k == "2")
				f2.load(file_name.c_str());
			if (*k == "13")
				f13.load(file_name.c_str());
			if (*k == "14")
				f14.load(file_name.c_str());
		} catch (JException& error) {
			FATAL(error.what() << endl);
		}

		NOTICE("OK" << endl);
	}

	double maximal_road_width = 0.0;  // [m]

	if (f1.intensity.rbegin()->getX() > maximal_road_width)
		maximal_road_width = f1.intensity.rbegin()->getX();
	if (f2.intensity.rbegin()->getX() > maximal_road_width)
		maximal_road_width = f2.intensity.rbegin()->getX();
	if (f13.intensity.rbegin()->getX() * getSinThetaC() > maximal_road_width)
		maximal_road_width = f13.intensity.rbegin()->getX() * getSinThetaC();
	if (f14.intensity.rbegin()->getX() * getSinThetaC() > maximal_road_width)
		maximal_road_width = f14.intensity.rbegin()->getX() * getSinThetaC();

	NOTICE("Maximal road width [m] " << maximal_road_width << endl);

	typedef JFunction1D_t::argument_type argument_type;
	typedef JFunction1D_t::result_type result_type;

	typedef JCDFTable1D<argument_type, result_type> JCDF1D_t;
	typedef JCDFTable2D<argument_type, result_type> JCDF2D_t;

	NOTICE("Setting up fast CDF tables... " << flush);

	JCDF1D_t f1s(f1, 200, 1.7);
	JCDF1D_t f2s(f2, 200, 1.7);
	JCDF2D_t f13s(f13, 200, 1.7);
	JCDF2D_t f14s(f14, 200, 1.7);

	struct {
		JCDF1D_t* integral;
		JCDF4D_t* function;
	} CDF[] = { { &f1s, &f1 }, { &f2s, &f2 } };

	struct {
		JCDF2D_t* integral;
		JCDF5D_t* function;
	} CDG[] = { { &f13s, &f13 }, { &f14s, &f14 } };

	NOTICE("OK" << endl);

	vector<JSharedPointer<JRadiationInterface> > radiation;

	if (true) {

		NOTICE("Setting up radiation tables... " << flush);

		JRadiation hydrogen(1.0, 1.0, 40, 0.01, 0.1, 0.1);
		JRadiation oxygen(8.0, 16.0, 40, 0.01, 0.1, 0.1);
		JRadiation chlorine(17.0, 35.0, 40, 0.01, 0.1, 0.1);

		JSharedPointer<JRadiation> Hydrogen(
				new JRadiationFunction(hydrogen, 300, 0.2, 1.0e11));
		JSharedPointer<JRadiation> Oxygen(
				new JRadiationFunction(oxygen, 300, 0.2, 1.0e11));
		JSharedPointer<JRadiation> Chlorine(
				new JRadiationFunction(chlorine, 300, 0.2, 1.0e11));

		JRadiationSource::source_type EErad(&JRadiation::TotalCrossSectionEErad,
				&JRadiation::EfromEErad);
		JRadiationSource::source_type Brems(&JRadiation::TotalCrossSectionBrems,
				&JRadiation::EfromBrems);
		JRadiationSource::source_type GNrad(&JRadiation::TotalCrossSectionGNrad,
				&JRadiation::EfromGNrad);

		radiation.push_back(
				new JRadiationSource(Oxygen, DENSITY_SEA_WATER * JSeaWater::O,
						EErad));
		radiation.push_back(
				new JRadiationSource(Chlorine,
						DENSITY_SEA_WATER * JSeaWater::Cl, EErad));
		radiation.push_back(
				new JRadiationSource(Hydrogen, DENSITY_SEA_WATER * JSeaWater::H,
						EErad));

		radiation.push_back(
				new JRadiationSource(Oxygen, DENSITY_SEA_WATER * JSeaWater::O,
						Brems));
		radiation.push_back(
				new JRadiationSource(Chlorine,
						DENSITY_SEA_WATER * JSeaWater::Cl, Brems));
		radiation.push_back(
				new JRadiationSource(Hydrogen, DENSITY_SEA_WATER * JSeaWater::H,
						Brems));

		radiation.push_back(
				new JRadiationSource(Oxygen, DENSITY_SEA_WATER * JSeaWater::O,
						GNrad));
		radiation.push_back(
				new JRadiationSource(Chlorine,
						DENSITY_SEA_WATER * JSeaWater::Cl, GNrad));
		radiation.push_back(
				new JRadiationSource(Hydrogen, DENSITY_SEA_WATER * JSeaWater::H,
						GNrad));

		NOTICE("OK" << endl);
	}

	JDetector detector;

	try {

		NOTICE("Load detector... " << flush);

		load(detectorFile, detector);

		NOTICE("OK" << endl);
	} catch (JException& error) {
		FATAL(error);
	}

	Header header;

	try {
		header = inputFile.getHeader();
	} catch (JException& error) {
		FATAL(error);
	}

	if (header.start_run) {
		Event::set_runNumber(header.start_run.run_id);
	}

	if (!header.coord_origin && header.can) {

		WARNING("Setting coordinate origin relative to bottom of can." << endl);

		header.coord_origin = MONTE_CARLO::coord_origin(0.0, 0.0,
				-header.can.zmin);
	}

	if (header.coord_origin) {

		detector -= JPosition3D(header.coord_origin.x, header.coord_origin.y,
				header.coord_origin.z);
		Zbed -= header.coord_origin.z;
	}

	JCylinder3D cylinder(detector.begin(), detector.end());

	cylinder.addMargin(maximal_road_width);

	if (cylinder.getZmin() < Zbed) {
		cylinder.setZmin(Zbed);
	}

	outputFile.open();

	if (!outputFile.is_open()) {
		FATAL("Error opening file " << outputFile << endl);
	}

	outputFile.put(header);
	outputFile.put(*gRandom);

	TProfile cpu("cpu", NULL, 14, 1.0, 7.0);
	TH1D job("job", NULL, 400, 0.5, 400.5);
	TProfile ems("ems", NULL, 14, 1.0, 7.0);

	JTimer timer;

	cout << "\nRunning with " << omp_get_max_threads() << " thread(s)...\n";

	int event_num = 0;

	vector<JEvent> eventList;

	double total_IO_time = 0;

	double read_time_start = omp_get_wtime();

	for (JMultipleFileScanner<Event>& in = inputFile; in.hasNext();) {

		++event_num;

		STATUS("event: " << setw(10) << in.getCounter() << '\r');
		DEBUG(endl);

		job.Fill(1.0);

		JEvent event(*in.next());

		event.removeHits();

		eventList.push_back(event);
	}

	double read_time_end = omp_get_wtime();

	total_IO_time += (read_time_end - read_time_start);

	double total_comp_time_start = omp_get_wtime();

#pragma omp parallel
	{
		TRandom *r3 = new TRandom3();
		r3->SetSeed(seed);

#pragma omp for schedule(dynamic, 1)
		for (int event_count = 0; event_count < 1000; ++event_count) {

			//timer.reset();
			//timer.start();

			for (vector<McTrack>::const_iterator track =
					eventList[event_count].TrackList().begin();
					track < eventList[event_count].TrackList().end(); ++track) {

				// (anti) muon

				if (is_muon(*track)) {

					job.Fill(2.0);

					pair<double, double> intersection;

					intersection = cylinder.getIntersection(getAxis(*track));

					double Zmin = intersection.first;
					double Zmax = intersection.second;

					if (Zmax <= 0.0) {
						continue;
					}

					// muon propagator

					JVertex vertex(0.0, track->t(), track->E());

					if (vertex.getZ() < Zbed) {

						// propagate muon through rock

						if (track->direction().z <= 0.0) {
							continue;
						}

						double ds = (Zbed - vertex.getZ())
								/ track->direction().z;

						if (vertex.getRange(gRock) > ds)
							vertex.step(gRock, ds);
						else
							continue;
					}

					if (vertex.getZ() < Zmin) {

						// propagate muon through water

						double ds = Zmin - vertex.getZ();

						if (vertex.getRange(gWater) > ds)
							vertex.step(gWater, ds);
						else
							continue;
					}

					if (vertex.getRange() <= 0.1) {
						continue;
					}

					job.Fill(3.0);

					JDetectorSubset subdetector(detector, getAxis(*track),
							maximal_road_width);

					if (subdetector.empty()) {
						continue;
					}

					job.Fill(4.0);

					JTrack muon(vertex);

					while (vertex.getE() >= Emin && vertex.getZ() < Zmax) {

						int N = radiation.size();

						double li[N];     // inverse interaction lengths
						double ls = 1.0e-5; // minimal total inverse interaction length (100 km)^-1

						for (int i = 0; i != N; ++i) {
							ls += li[i] =
									radiation[i]->getInverseInteractionLength(
											vertex.getE());
						}

						double temp;

						temp = r3->Exp(1.0);

						double ds = min(temp / ls, vertex.getRange());

						vertex.step(ds);

						if (vertex.getE() >= Emin) {

							double Es = 0.0;      // shower energy [GeV]
							double y;

							y = r3->Uniform(ls);

							for (int i = 0; i != N; ++i) {

								y -= li[i];

								if (y < 0.0) {
									Es = radiation[i]->getEnergyOfShower(
											vertex.getE()); // shower energy [GeV]
									break;
								}
							}

							vertex.applyEloss(Es);

							if (Es >= Ecut) {
								muon.push_back(vertex);
							}
						}
					}

					if (vertex.getE() < Emin && vertex.getRange() > 0.0) {

						vertex.step(vertex.getRange());

						muon.push_back(vertex);
					}

					//double hit_time_start = omp_get_wtime();

					// generate direct and scattered light from muon

					Zmin = muon.begin()->getZ()
							- maximal_road_width / getTanThetaC();
					Zmax = muon.rbegin()->getZ();

					for (JDetector::const_iterator module =
							subdetector.lower_bound(Zmin);
							module < subdetector.lower_bound(Zmax); ++module) {

						double z0 = muon.begin()->getZ();
						double t0 = muon.begin()->getT();
						double R = module->getX();

						if (module->getZ() - R / getTanThetaC()
								>= muon.begin()->getZ()
								&& module->getZ() - R / getTanThetaC()
										<= muon.rbegin()->getZ()) {

							for (int i = 0; i != 2; ++i) {

								// for (i = 0)
								double Rmax = CDF[i].integral->rbegin()->getX();

								if (R < Rmax) {

									try {

										double NPE;

										NPE = CDF[i].integral->getNPE(R)
												* module->size();

										int N;

										N = r3->Poisson(NPE);

										double ns = 0.0;

										for (JModule::const_iterator pmt =
												module->begin();
												pmt != module->end() && N != 0;
												++pmt) {

											double R = pmt->getX();
											double Z = z0 - pmt->getZ();
											double theta = pmt->getTheta();
											double phi = fabs(pmt->getPhi());

											double npe;

#pragma omp critical
											{
												npe = CDF[i].function->getNPE(R,
														theta, phi);
											}

											ns += npe;

											int n1 = getNumberOfPhotoElectrons(
													NPE, N, npe, &(*r3));

											job.Fill((double) (101 + i),
													(double) n1);

											while (n1 != 0) {

												double temp;

												temp = r3->Rndm();

												double dt;

#pragma omp critical
												{
													dt =
															CDF[i].function->getTime(
																	R, theta,
																	phi, temp);
												}

												int n =
														getNumberOfPhotoElectrons(
																n1);

												double t1 = t0
														+ (R * getTanThetaC()
																- Z) / C + dt;
												double a1 = n;
												int type =
														(i == 0 ?
																HIT_TYPE_MUON_DIRECT :
																HIT_TYPE_MUON_SCATTERED);
												uint id =
														eventList[event_count].HitList().size()
																+ 1;

												McHit hit(id, pmt->getID(), t1,
														a1, type, track->id(),
														a1, t1);

												eventList[event_count].Add(hit);

												n1 -= n;
											}
										}

										if (ns > NPE) {

											job.Fill((double) (301 + i));

										}
									} catch (exception& error) {

										job.Fill((double) (201 + i));

									}
								}
							}
						}
					}

					// generate direct and scattered light from EM showers.

					for (JTrack::const_iterator vertex = muon.begin();
							vertex < muon.end(); ++vertex) {

						double t0 = vertex->getT();
						double z0 = vertex->getZ();
						double Es = vertex->getEs();

						if (Es >= Ecut) {

							int trackID = track->id();

							if (writeEMShowers) {
								trackID =
										eventList[event_count].secondary().size()
												+ 1;
							}

							bool secondaryHit = false;

							Zmin = z0 - maximal_road_width;
							Zmax = z0 + maximal_road_width;

							for (JDetector::const_iterator module =
									subdetector.lower_bound(Zmin);
									module < subdetector.lower_bound(Zmax);
									++module) {

								double R = module->getX();
								double Z = z0 - module->getZ();

								double D = sqrt(R * R + Z * Z);
								double cd = -Z / D;

								for (int i = 0; i != 2; ++i) {

									double Dmax =
											CDG[i].integral->rbegin()->getX();

									if (D < Dmax) {

										try {

											double NPE;

											NPE = CDG[i].integral->getNPE(D, cd)
													* Es * module->size();

											int N;

											N = r3->Poisson(NPE);

											double ns = 0.0;

											for (JModule::const_iterator pmt =
													module->begin();
													pmt < module->end();
													++pmt) {

												if (N != 0) {

													double R = pmt->getX();
													double Z = z0 - pmt->getZ();
													double theta =
															pmt->getTheta();
													double phi = fabs(
															pmt->getPhi());
													double D = sqrt(
															R * R + Z * Z);
													double cd = -Z / D;

													double npe;

#pragma omp critical
													{
														npe =
																CDG[i].function->getNPE(
																		D, cd,
																		theta,
																		phi)
																		* Es;
													}

													ns += npe;

													int n1 =
															getNumberOfPhotoElectrons(
																	NPE, N, npe,
																	&(*r3));

													job.Fill((double) (113 + i),
															(double) n1);

													while (n1 != 0) {

														double temp;

														temp = r3->Rndm();

														double z =
																Z
																		+ geanz.getLength(
																				Es,
																				temp);
														double D = sqrt(
																R * R + z * z);
														double cd = -z / D;

														temp = r3->Rndm();

														double dt;

#pragma omp critical
														{
															dt =
																	CDG[i].function->getTime(
																			D,
																			cd,
																			theta,
																			phi,
																			temp);
														}

														int n =
																getNumberOfPhotoElectrons(
																		n1);

														double t1 =
																t0
																		+ (z - Z
																				+ D
																						* getIndexOfRefraction())
																				/ C
																		+ dt;
														double a1 = n;
														int type =
																(i == 0 ?
																		HIT_TYPE_SHOWER_DIRECT :
																		HIT_TYPE_SHOWER_SCATTERED);
														uint id =
																eventList[event_count].HitList().size()
																		+ 1;

														McHit hit(id,
																pmt->getID(),
																t1, a1, type,
																trackID, a1,
																t1);

														eventList[event_count].Add(
																hit);

														n1 -= n;

														secondaryHit = true;
													}
												}
											}

											if (ns > NPE) {

												job.Fill((double) (313 + i));

											}
										} catch (exception& error) {

											job.Fill((double) (213 + i));

										}
									}
								}
							}

							if (writeEMShowers && secondaryHit) {

								JPosition3D pos(0.0, 0.0, z0);
								JDirection3D dir = getDirection(*track);

								pos.add(getPosition(*track));
								pos.rotate_back(JRotation3D(dir));

								Secondary showerOnTrack;

								showerOnTrack.set_id(trackID);
								showerOnTrack.set_t(t0);
								showerOnTrack.set_position(pos.getX(),
										pos.getY(), pos.getZ());
								showerOnTrack.set_direction(dir.getDX(),
										dir.getDY(), dir.getDZ());
								showerOnTrack.set_E(Es);
								showerOnTrack.set_type(-1); // negative gamma particle code
								showerOnTrack.set_origin(track->id());
								showerOnTrack.validate();

								eventList[event_count].Add(showerOnTrack);

							}
						}
					}

					//double hit_time_end = omp_get_wtime();

					//total_hit_time += (hit_time_end - hit_time_start);

					Zmin = muon.begin()->getZ();
					Zmax = muon.rbegin()->getZ();

					if (Zmax > Zmin) {

						ems.Fill(log10(muon.begin()->getE()),
								(double) muon.size() / (Zmax - Zmin));

					}
				}
			}

			// merge hits on same PMT within maximal time window

			eventList[event_count].mergeHits(Tmax);

			// timer.stop();

			if (eventList[event_count].hasNeutrino()
					&& !eventList[event_count].HitList().empty()) {

				cpu.Fill(log10(eventList[event_count].neutrino().E()),
						(double) timer.usec_ucpu * 1.0e-3);

			}
		}
	}

	double total_comp_time_end = omp_get_wtime();

	double write_time_start = omp_get_wtime();

	for (int event_count = 0; event_count < 1000; ++event_count) {

		if (!eventList[event_count].HitList().empty()) {

			outputFile.put(eventList[event_count]);

			job.Fill(10.0);

		}

	}

	double write_time_end = omp_get_wtime();

	total_IO_time += (write_time_end - write_time_start);

	STATUS(endl);

	outputFile.put(*gRandom);

	outputFile.close();

	/*cout << "Total time - hit detection - all events = " << total_hit_time
	 << " seconds.\n";*/
	cout << "Total time - all events = "
			<< total_comp_time_end - total_comp_time_start << " seconds.\n";
	cout << "Total time - IO = " << total_IO_time << " seconds.\n";
}
