
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
#include "JPhysics/JPDFLibrary.hh"
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


/**   --------------------------------------------------------------------------
 *   |                                                                          |
 *   |  Program to simulate detector response to muons and showers.             |
 *   |                                                                          |
 *   |   - PDF file descriptor should contain '%' character;                    |
 *   |     The file names are obtained by replacing % with 1, 2, 13 and 14;     |
 *   |                                                                          |
 *   |   - Detector geometry is read from .det file;                            |
 *   |                                                                          |
 *   |   - Events are read from / written to ROOT formatted .evt file;          |
 *    --------------------------------------------------------------------------
 */
int main(int argc, char **argv)
{
  using namespace std;

  string         fileDescriptor;
  vector<string> inputFile;
  string         outputFile;
  string         interfaceInput;
  Long64_t       numberOfEvents;
  double         Tmax;
  bool           geasim;
  bool           writeEMShowers;
  int            debug;

  try { 

    JParser<> zap;
    
    zap['F'] = make_field(fileDescriptor);
    zap['f'] = make_field(inputFile);
    zap['o'] = make_field(outputFile)        = "sirene.root";
    zap['a'] = make_field(interfaceInput);
    zap['n'] = make_field(numberOfEvents)    = numeric_limits<Long64_t>::max();
    zap['T'] = make_field(Tmax)              = 0.0;
    zap['G'] = make_field(geasim);
    zap['s'] = make_field(writeEMShowers);
    zap['d'] = make_field(debug)             = 1;
    
    if (zap.read(argc, argv) != 0)
      return 1;
  }
  catch(const exception &error) {
    ERROR(error.what() << endl);
    return 2;
  }


  using namespace JTOOLS;
  using namespace JPHYSICS;
  using namespace JSIRENE;
  using namespace JEEP;
  using namespace JLANG;

  typedef JSplineFunction1D_t                                     JFunction1D_t;
  typedef 
    JMapList<JPolint1FunctionalMap,
    JMapList<JPolint1FunctionalGridMap,
    JMapList<JPolint1FunctionalGridMap> > >                       J3DMap_t;
  typedef 
    JMapList<JPolint1FunctionalMap,
    JMapList<JPolint1FunctionalMap,
    JMapList<JPolint1FunctionalGridMap,
    JMapList<JPolint1FunctionalGridMap> > > >                     J4DMap_t;

  typedef JCDFTable<JFunction1D_t, J3DMap_t>                      JCDF4D_t;
  typedef JCDFTable<JFunction1D_t, J4DMap_t>                      JCDF5D_t;

  JCDF4D_t f1;       // direct    light from muon
  JCDF4D_t f2;       // scattered light from muon
  JCDF5D_t f13;      // direct    light from EM shower
  JCDF5D_t f14;      // scattered light from EM shower

  const double a    = 2.67e-1 * DENSITY_SEA_WATER;     // Ionisation energy loss [Gev/m]
  const double Ecut = 0.1;                             // minimal energy for generation of light from shower [GeV]
  const double Emin = 1.0;                             // minimal energy of muon for shower generation [GeV]
  double       Zbed = 0.0;                             // level of seabed [m]


  for (string buffer[] = { "1", "2", "13", "14", ""}, *k = buffer; *k != ""; ++k) {

    string file_name = fileDescriptor;

    string::size_type ipos = file_name.find('%');

    if (ipos == string::npos) FATAL("error file descriptor " << fileDescriptor << endl);

    file_name.replace(ipos, 1, *k);
    
    NOTICE("loading input from file " << file_name << "... " << flush);

    try {
      if (*k ==  "1") f1 .load(file_name.c_str());
      if (*k ==  "2") f2 .load(file_name.c_str());
      if (*k == "13") f13.load(file_name.c_str());
      if (*k == "14") f14.load(file_name.c_str());
    }      
    catch(const JException& error) {
      FATAL(error.what() << endl);
    }
    
    NOTICE("OK" << endl);
  }


  double maximal_road_width = 0.0;  // [m]

  if (f1 .intensity.rbegin()->first > maximal_road_width) maximal_road_width = f1 .intensity.rbegin()->first;
  if (f2 .intensity.rbegin()->first > maximal_road_width) maximal_road_width = f2 .intensity.rbegin()->first;
  if (f13.intensity.rbegin()->first > maximal_road_width) maximal_road_width = f13.intensity.rbegin()->first;
  if (f14.intensity.rbegin()->first > maximal_road_width) maximal_road_width = f14.intensity.rbegin()->first;

  NOTICE("Maximal road width [m] " << maximal_road_width << endl);

  NOTICE("Setting up fast CDF tables... " << flush);

  typedef JFunction1D_t::argument_type               argument_type;
  typedef JFunction1D_t::result_type                 result_type;

  typedef JCDFTable1D<argument_type, result_type>    JCDF1D_t;
  typedef JCDFTable2D<argument_type, result_type>    JCDF2D_t;
  
  const JCDF1D_t f1s (f1,  200);
  const JCDF1D_t f2s (f2,  200);
  const JCDF2D_t f13s(f13, 200);
  const JCDF2D_t f14s(f14, 200);

  struct {
    const JCDF1D_t* integral;
    const JCDF4D_t* function;
  } CDF[] = { { &f1s, &f1 }, { &f2s, &f2 } };

  struct {
    const JCDF2D_t* integral;
    const JCDF5D_t* function;
  } CDG[] = { { &f13s, &f13 }, { &f14s, &f14 } };

  NOTICE("OK" << endl);


  vector< JSharedPointer<JRadiationInterface> > radiation;

  if (true) {

    NOTICE("Setting up radiation tables... " << flush);

    const JRadiation hydrogen( 1.0,  1.0, 40, 0.01, 0.1, 0.1);
    const JRadiation oxygen  ( 8.0, 16.0, 40, 0.01, 0.1, 0.1);
    const JRadiation chlorine(17.0, 35.0, 40, 0.01, 0.1, 0.1);

    const JRadiationFunction Hydrogen(hydrogen, 300, 0.2, 1.0e11);
    const JRadiationFunction Oxygen  (oxygen,   300, 0.2, 1.0e11);
    const JRadiationFunction Chlorine(chlorine, 300, 0.2, 1.0e11);
    
    JRadiationSource::source_type EErad(&JRadiation::TotalCrossSectionEErad, &JRadiation::EfromEErad);
    JRadiationSource::source_type Brems(&JRadiation::TotalCrossSectionBrems, &JRadiation::EfromBrems);
    JRadiationSource::source_type GNrad(&JRadiation::TotalCrossSectionGNrad, &JRadiation::EfromGNrad);

    radiation.push_back(JRadiationSource(Oxygen,   DENSITY_SEA_WATER * JSeaWater::O,  EErad));
    radiation.push_back(JRadiationSource(Chlorine, DENSITY_SEA_WATER * JSeaWater::Cl, EErad));
    radiation.push_back(JRadiationSource(Hydrogen, DENSITY_SEA_WATER * JSeaWater::H,  EErad));

    radiation.push_back(JRadiationSource(Oxygen,   DENSITY_SEA_WATER * JSeaWater::O,  Brems));
    radiation.push_back(JRadiationSource(Chlorine, DENSITY_SEA_WATER * JSeaWater::Cl, Brems));
    radiation.push_back(JRadiationSource(Hydrogen, DENSITY_SEA_WATER * JSeaWater::H,  Brems));

    radiation.push_back(JRadiationSource(Oxygen,   DENSITY_SEA_WATER * JSeaWater::O,  GNrad));
    radiation.push_back(JRadiationSource(Chlorine, DENSITY_SEA_WATER * JSeaWater::Cl, GNrad));
    radiation.push_back(JRadiationSource(Hydrogen, DENSITY_SEA_WATER * JSeaWater::H,  GNrad));

    NOTICE("OK" << endl);
  }

  JDetector detector;

  try {
    load(interfaceInput, detector);
  }
  catch(const JException& error) {
    FATAL(error);
  }

  // ROOT I/O

  for (vector<string>::const_iterator i = inputFile.begin(); i != inputFile.end(); ++i)
    MonteCarloEvent_Reader.Add(i->c_str());

  Header header(MonteCarloEvent_Reader);

  NOTICE("Monte Carlo Event tree: " << MonteCarloEvent_Reader.GetEntries() << " " << header.getNumberOfHeaders() << endl);
  
  if (!header.can) {

    WARNING("Recovering can." << endl);

    JCylinder cyl(detector);

    header.can = MONTE_CARLO::can(cyl.getZmin(), cyl.getZmax(), cyl.getRadius());
    
    // expand can

    header.can.zmin -= maximal_road_width;
    header.can.zmax += maximal_road_width;
    header.can.r    += maximal_road_width;
  }

  if (header.start_run) Event::set_runNumber(header.start_run.run_id);

  if (!header.coord_origin && header.can) {

    WARNING("Setting coordinate origin relative to bottom of can." << endl);

    header.coord_origin = MONTE_CARLO::coord_origin(0.0, 0.0, -header.can.zmin);
  }

  if (header.coord_origin) {

    detector -= JPoint3D(header.coord_origin.x,
			 header.coord_origin.y,
			 header.coord_origin.z);
    Zbed     -=          header.coord_origin.z;

  } else
    WARNING("Missing coordinate origin." << endl);


  TFile* out = new TFile(outputFile.c_str(), "recreate");

  if (out == NULL || !out->IsOpen()) FATAL("Error opening file " << outputFile << endl);


  header.Write();
  
  MonteCarloEvent_Writer.SetDirectory(out);

  TProfile cpu("cpu", NULL,  14,  1.0,   7.0);
  TH1D     job("job", NULL, 300,  0.5, 300.5);

  JTimer timer;


  numberOfEvents = min(numberOfEvents, MonteCarloEvent_Reader.GetEntries());

  for (Long64_t event_count = 0; event_count < numberOfEvents; ++event_count) {

    STATUS("event: " << setw(10) << event_count << '\r'); DEBUG(endl);

    MonteCarloEvent_Reader.GetEvent(event_count);

    job.Fill(1.0);

    JEvent event(*MonteCarloEvent_Reader.GetAddress());

    event.removeHits();

    timer.reset();
    timer.start();

    for (vector<McTrack>::const_iterator track = event.TrackList().begin(); track != event.TrackList().end(); ++track) {

      // (anti) muon

      if (is_muon(*track)) {

	job.Fill(2.0);

	double E0 = track->E();
	double z0 = 0.0;
	double t0 = track->t();
      
	const pair<double, double> intersection = getIntersection(*track, header.can);

	const double Zmin = intersection.first;
	const double Zmax = intersection.second;

	if (Zmax <= 0.0)       continue;

	if (Zbed > track->position().z) {

	  // propagate muon through rock
	    
	  const double ds = (Zbed - track->position().z) / track->direction().z;

	  if (geaneRock(E0) < ds)  continue;

	  E0  = geaneRock(E0, ds);
	  z0 += ds;
	  t0 += ds / C;
	}

	if (Zmin > z0) {

	  // propagate muon through water
	    
	  const double ds = Zmin - z0;
	  
	  if (geaneWater(E0) < ds)  continue;

	  E0  = geaneWater(E0, ds);
	  z0 += ds;
	  t0 += ds / C;
	}

	job.Fill(3.0);

	const JGeometry trackGeometry = getGeometry(*track);
        const JRotation3D trackRotation(trackGeometry.getDirection());

        JPoint3D trackPosition = trackGeometry.getPosition();
        trackPosition.rotate(trackRotation);
        JDirection3D trackDirection = trackGeometry.getDirection();

	const JDetectorSubset subdetector(detector, getGeometry(*track), maximal_road_width);

	if (subdetector.empty()) continue;

	job.Fill(4.0);


	while (E0 > MASS_MUON * getIndexOfRefraction() && z0 < Zmax) {

	
	  const int N = radiation.size();
	
	  double li[N];                                       // inverse interaction lengths
	  double ls = 1.0e-5;                                 // minimal total inverse interaction length (100 km)^-1
	
	  for (int i = 0; i != N; ++i)
	    ls += li[i] = radiation[i]->getInverseInteractionLength(E0);

	  const double step = min(gRandom->Exp(1.0) / ls, E0/a);

      
	  E0 -= a*step;                                       // ionisation energy loss
	
	  double Es = 0.0;                                    // shower energy [GeV]

	  if (E0 >= Emin) {

	    double y  = gRandom->Uniform(ls);

	    for (int i = 0; i != N; ++i) {
	      
	      y -= li[i];
	      
	      if (y < 0.0) {
		Es = radiation[i]->getEnergyOfShower(E0);     // shower energy [GeV]
		break;
	      }
	    }
	  }


	  // generate direct and scattered light from muon
	  
	  for (JDetector::const_iterator module = subdetector.lower_bound(z0); module != subdetector.end(); ++module) {

	    const double R = module->getX();
	    const double Z = z0 - module->getZ();

	    if (Z + step <  -maximal_road_width / getTanThetaC()) break;

	    if (Z        <= -R / getTanThetaC() &&
		Z + step >  -R / getTanThetaC()) {

	      for (int i = 0; i != 2; ++i) {
		      
		const double Rmax = CDF[i].integral->getUpperKey();
		
		if (R < Rmax) {
		
		  try {

		    const double NPE = CDF[i].integral->getNPE(R) * module->size();
		    const int    N   = gRandom->Poisson(NPE); 
		    
		    if (N != 0) {

		      for (JModule::const_iterator pmt = module->begin(); pmt != module->end(); ++pmt) {

			const double R     = pmt->getX();
			const double Z     = z0 - pmt->getZ();
			const double theta = pmt->getTheta();
			const double phi   = pmt->getPhi();
			
			const double npe = CDF[i].function->getNPE(R, theta, phi);

			int n1 = getNumberOfPhotoElectrons(NPE, N, npe);

			job.Fill((double) (101 + i), (double) n1);

			while (n1 != 0) {
			  
			  const double dt = CDF[i].function->getTime(R, theta, phi, gRandom->Rndm());
			  const int    n  = getNumberOfPhotoElectrons(n1);		  
			  
			  const double t1   =  t0  +  (R * getTanThetaC() - Z) / C  +  dt; 
			  const double a1   =  n;
			  const int    type = (i == 0 ? +5 : -5);
			  const uint   id   =  event.HitList().size() + 1;
			  
			  const Hit hit(id, pmt->getID(), t1, a1, type, track->id(), a1, t1);
			  
			  event.Add(hit);

			  n1 -= n;
			}
		      }
		    }
		  }
		  catch(const exception& error) {
		    job.Fill((double) (201 + i));
		  }
		}
	      }
	    }
	  }
	  

	  // generate direct and scattered light from EM shower.

	  if (Es >= Ecut) {
	    
	    // store shower at z0+step as Secondary with energy Es and code -1
	    if(writeEMShowers) {
	      JPoint3D showerPosition(trackPosition.getX(),trackPosition.getY(),z0+step);
	      showerPosition.rotate_back(trackRotation);
	      
	      Secondary showerOnTrack;
	      showerOnTrack.set_id(event.secondary().size()+1);
	      showerOnTrack.set_t(t0+step/C);
	      showerOnTrack.set_position(showerPosition.getX(),showerPosition.getY(),showerPosition.getZ());
	      showerOnTrack.set_direction(trackDirection.getDX(),trackDirection.getDY(),trackDirection.getDZ());
	      showerOnTrack.set_E(Es);
	      showerOnTrack.set_type(-1); //negative gamma particle code
	      showerOnTrack.set_origin(track->id());
	      showerOnTrack.validate();
	      
	      event.Add(showerOnTrack);
	    }
	    
	    for (JDetector::const_iterator module = subdetector.lower_bound(z0 - maximal_road_width); module != subdetector.end(); ++module) {
	      
	      const double R = module->getX();
	      const double Z = z0 - module->getZ();
	      
	      if (Z < -maximal_road_width) break;

	      const double D  = sqrt(R*R + Z*Z);
	      const double cd = -Z / D; 

	      for (int i = 0; i != 2; ++i) {
	    
		const double Dmax = CDG[i].integral->getUpperKey();
	    
		if (D < Dmax) {
		
		  try {

		    const double NPE = CDG[i].integral->getNPE(D, cd) * Es * module->size();
		    const int    N   = gRandom->Poisson(NPE); 

		    if (N != 0) {

		      for (JModule::const_iterator pmt = module->begin(); pmt != module->end(); ++pmt) {

			const double R     = pmt->getX();
			const double Z     = z0 - pmt->getZ();
			const double D     = sqrt(R*R + Z*Z);
			const double cd    = -Z / D; 
			const double theta = pmt->getTheta();
			const double phi   = pmt->getPhi();
			    
			const double npe = CDG[i].function->getNPE(D, cd, theta, phi) * Es;
		  
			int n1 = getNumberOfPhotoElectrons(NPE, N, npe);

			job.Fill((double) (113 + i), (double) n1);

			while (n1 != 0) {

			  const double z  =  Z  +  geanz.getLength(Es, gRandom->Rndm());
			  const double D  =  sqrt(R*R + z*z);
			  const double cd = -z / D; 
			  
			  const double dt = CDG[i].function->getTime(D, cd, theta, phi, gRandom->Rndm());
			  const int    n  = getNumberOfPhotoElectrons(n1);		  
			  
			  const double t1   = t0  +  (z - Z +  D * getIndexOfRefraction()) / C  +  dt;
			  const double a1   = n;
			  const int    type = i == 0 ? +3 : -3;
			  const uint   id   = event.HitList().size() + 1;
			  
			  const Hit hit(id, pmt->getID(), t1, a1, type, track->id(), a1, t1);
			  
			  event.Add(hit);
			  
			  n1 -= n;
			}
		      }
		    }
		  }
		  catch(const exception& error) {
		    job.Fill((double) (213 + i));
		  }
		}
	      }
	    }
	  }

	  E0 -= Es;
	  z0 += step;
	  t0 += step / C;
	}
      }
    }


    if (geasim) {

      if (event.hasNeutrino()) {

	// One particle approximation for all secondaries from the neutrino interaction vertex.

	const Neutrino& neutrino = event.neutrino();
	const Vec3D&    vertex   = neutrino.position();
	const double    R        = sqrt(vertex.x*vertex.x + vertex.y*vertex.y);

	if (vertex.z >= header.can.zmin && 
	    vertex.z <= header.can.zmax && 
	    R        <= header.can.r) {

	  Vec3D shower;
	  
	  const double   z0 = 0.0;
	  const double   t0 = neutrino.t();
	  double         Es = 0.0;             // total  energy   [GeV]
	  vector<double> Ei;                   // shower energies [GeV]
      
	  for (vector<McTrack>::const_iterator track = event.TrackList().begin(); track != event.TrackList().end(); ++track) {
	  
	    if (!is_muon(*track)) {
	  
	      if (fabs(vertex.x - track->position().x) < 0.1 &&
		  fabs(vertex.y - track->position().y) < 0.1 &&
		  fabs(vertex.z - track->position().z) < 0.1) {
	      
		const double E = track->E() * phytia(track->type());

		Es     += E;
		shower += E * track->direction();

		Ei.push_back(E);
	      }
	    }
	  }
      
	  if (Es >= Ecut) {
	
	    // generate direct and scattered light from secondary particles.
	  
	    const JDetectorSubset subdetector(detector, JGeometry(getPosition(vertex), getDirection(shower)), maximal_road_width);
	  
	    for (JDetector::const_iterator module = subdetector.lower_bound(z0 - maximal_road_width); module != subdetector.end(); ++module) {

	      const double R  = module->getX();
	      const double Z  = z0 - module->getZ();

	      if (Z < -maximal_road_width) break;

	      const double D  = sqrt(R*R + Z*Z);
	      const double cd = -Z / D; 
	      
	      for (int i = 0; i != 2; ++i) {
		
		const double Dmax = CDG[i].integral->getUpperKey();
		
		if (D < Dmax) {
		  
		  try {
		    
		    const double NPE = CDG[i].integral->getNPE(D, cd) * Es * module->size();
		    const int    N   = gRandom->Poisson(NPE); 
		    
		    if (N != 0) {
		      
		      for (JModule::const_iterator pmt = module->begin(); pmt != module->end(); ++pmt) {
			
			const double R     = pmt->getX();
			const double Z     = z0 - pmt->getZ();
			const double D     = sqrt(R*R + Z*Z);
			const double cd    = -Z / D; 
			const double theta = pmt->getTheta();
			const double phi   = pmt->getPhi();
			
			const double npe = CDG[i].function->getNPE(D, cd, theta, phi) * Es;
			
			int n1 = getNumberOfPhotoElectrons(NPE, N, npe);

			job.Fill((double) (123 + i), (double) n1);

			while (n1 != 0) {
			  
			  double z  = 0.0;
			  double y  = gRandom->Uniform(Es);

			  for (vector<double>::const_iterator energy = Ei.begin(); energy != Ei.end(); ++energy) {
			    
			    y -= *energy;
			    
			    if (y < 0.0) {
			      z = Z  +  geanz.getLength(*energy, gRandom->Rndm());
			      break;
			    }
			  }

			  const double D  =  sqrt(R*R + z*z);
			  const double cd = -z / D; 
			
			  const double dt = CDG[i].function->getTime(D, cd, theta, phi, gRandom->Rndm());
			  const int    n  = getNumberOfPhotoElectrons(n1);		  
			
			  const double t1   = t0  +  (z - Z +  D * getIndexOfRefraction()) / C  +  dt;
			  const double a1   = n;
			  const int    type = i == 0 ? +99 : -99;
			  const uint   id   = event.HitList().size() + 1;
			
			  const Hit hit(id, pmt->getID(), t1, a1, type, 0, a1, t1);
			
			  event.Add(hit);
			
			  n1 -= n;
			}
		      }
		    }
		  }
		  catch(const exception& error) {
		    job.Fill((double) (223 + i));
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    // merge hits on same PMT within maximal time window

    event.mergeHits(Tmax);

    timer.stop();

    if (event.hasNeutrino() && !event.HitList().empty())
      cpu.Fill(log10(event.neutrino().E()), (double) timer.usec_ucpu * 1.0e-3);

    if (!event.HitList().empty()) {

      MonteCarloEvent_Writer.Write(event);

      job.Fill(5.0);
    }
  }
  STATUS(endl); 

  NOTICE("Job summary" << endl);

  NOTICE(left << setw(30) << "Number of events input      " << right << setw(10) << (int) job.GetBinContent(  1) << endl);
  NOTICE(left << setw(30) << "Number of muons             " << right << setw(10) << (int) job.GetBinContent(  2) << endl);
  NOTICE(left << setw(30) << "Number of muons in can      " << right << setw(10) << (int) job.GetBinContent(  3) << endl);
  NOTICE(left << setw(30) << "Number of muons within road " << right << setw(10) << (int) job.GetBinContent(  4) << endl);
  NOTICE(left << setw(30) << "Number of events output     " << right << setw(10) << (int) job.GetBinContent(  5) << endl);

  NOTICE(left << setw(30) << "Number of photons  1        " << right << setw(10) << (int) job.GetBinContent(101) << endl);
  NOTICE(left << setw(30) << "Number of photons  2        " << right << setw(10) << (int) job.GetBinContent(102) << endl);
  NOTICE(left << setw(30) << "Number of photons 13        " << right << setw(10) << (int) job.GetBinContent(113) << endl);
  NOTICE(left << setw(30) << "Number of photons 14        " << right << setw(10) << (int) job.GetBinContent(114) << endl);
  NOTICE(left << setw(30) << "Number of photons 23        " << right << setw(10) << (int) job.GetBinContent(123) << endl);
  NOTICE(left << setw(30) << "Number of photons 24        " << right << setw(10) << (int) job.GetBinContent(124) << endl);

  NOTICE(left << setw(30) << "Number of errors   1        " << right << setw(10) << (int) job.GetBinContent(201) << endl);
  NOTICE(left << setw(30) << "Number of errors   2        " << right << setw(10) << (int) job.GetBinContent(202) << endl);
  NOTICE(left << setw(30) << "Number of errors  13        " << right << setw(10) << (int) job.GetBinContent(213) << endl);
  NOTICE(left << setw(30) << "Number of errors  14        " << right << setw(10) << (int) job.GetBinContent(214) << endl);
  NOTICE(left << setw(30) << "Number of errors  23        " << right << setw(10) << (int) job.GetBinContent(223) << endl);
  NOTICE(left << setw(30) << "Number of errors  24        " << right << setw(10) << (int) job.GetBinContent(224) << endl);

  out->Write();
  out->Close();
}
