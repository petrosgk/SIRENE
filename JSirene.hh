#ifndef __JSIRENE__
#define __JSIRENE__

#include <vector>
#include <iterator>

#include "TRandom3.h"

#include "antcc/Header.hh"
#include "antcc/Track.hh" 
#include "antcc/Event.hh" 

#include "JTools/JConstants.hh"
#include "JPhysics/JPDFToolkit.hh"
#include "JGeometry3D/JPosition3D.hh"
#include "JGeometry3D/JVersor3D.hh"
#include "JGeometry3D/JAxis3D.hh"
#include "JGeometry3D/JTrack3D.hh"
#include "JDetector/JDetector.hh"
#include "JSirene/JAntcc.hh"
#include "JLang/JException.hh"

namespace JSIRENE {

//typedef  unsigned char    uchar;
//typedef  unsigned short   ushort;
//typedef  unsigned int     uint;

using JGEOMETRY3D::JPosition3D;
using JGEOMETRY3D::JVersor3D;
using JGEOMETRY3D::JAxis3D;
using JGEOMETRY3D::JTrack3D;
using JTOOLS::getSpeedOfLight;
using JTOOLS::MASS_MUON;
using JTOOLS::getIndexOfRefraction;
using JTOOLS::DENSITY_SEA_WATER;
using JPHYSICS::JGeane;
using JPHYSICS::gWater;
using JPHYSICS::geanz;
using JLANG::JEmptyCollection;

/**
 * Get position.
 *
 * \param  v        vector
 * \return          position
 */
inline JPosition3D getPosition(const Vec3D& v) {
	return JPosition3D(v.x, v.y, v.z);
}

/**
 * Get position.
 *
 * \param  track    track
 * \return          position
 */
inline JPosition3D getPosition(const BaseTrack& track) {
	return getPosition(track.position());
}

/**
 * Get direction.
 *
 * \param  v        vector
 * \return          direction
 */
inline JVersor3D getDirection(const Vec3D& v) {
	return JVersor3D(v.x, v.y, v.z);
}

/**
 * Get direction.
 *
 * \param  track    track
 * \return          direction
 */
inline JVersor3D getDirection(const BaseTrack& track) {
	return getDirection(track.direction());
}

/**
 * Get axis.
 *
 * \param  track    track
 * \return          axis
 */
inline JAxis3D getAxis(const BaseTrack& track) {
	return JAxis3D(getPosition(track), getDirection(track));
}

/**
 * Get track.
 *
 * \param  track    track
 * \return          track
 */
inline JTrack3D getTrack(const BaseTrack& track) {
	return JTrack3D(getAxis(track), track.t());
}

// GEANT particle codes

inline bool is_photon(const McTrack& track) {
	return track.type() == 1 || track.type() == 7;
}
inline bool is_electron(const McTrack& track) {
	return track.type() == 2 || track.type() == 3;
}
inline bool is_neutrino(const McTrack& track) {
	return track.type() == 4;
}
inline bool is_muon(const McTrack& track) {
	return track.type() == 5 || track.type() == 6;
}
inline bool is_pion(const McTrack& track) {
	return track.type() == 8 || track.type() == 9;
}
inline bool is_neutron(const McTrack& track) {
	return track.type() == 13;
}
inline bool is_proton(const McTrack& track) {
	return track.type() == 14;
}
inline bool is_tau(const McTrack& track) {
	return track.type() == 33 || track.type() == 34;
}

/**
 * Get number of photo-electrons of a hit given the expectation
 * values of the number of photo-electrons on a module and PMT.
 *
 * The return value is evaluated by pick-and-drop statistics from
 * the generated number of photo-electrons when the expectation
 * value of the number of photo-electrons on a module deviates
 * less than 5 sigmas from 0 (i.e. when it is less than 25).
 * Otherwise, the return value is evaluated by Poisson statistics
 * from the expectation value of the number of photo-electrons on PMT.
 *
 * \param  NPE    expectation value of npe on module
 * \param  N      generated   value of npe on module
 * \param  npe    expectation value of npe on PMT
 * \return        number of photo-electrons on PMT
 */
inline int getNumberOfPhotoElectrons(const double NPE, const int N,
		const double npe, TRandom *r3) {

	double temp;

	if (NPE < 25.0) {

		int n = 0;

		for (int j = N; j != 0; --j) {

			temp = r3->Rndm();

			if (temp * NPE < npe)
				++n;
		}

		return n;

	} else {

		temp = r3->Poisson(npe);

		return temp;
	}
}

/**
 * Get number of photo-electrons of a hit given number of photo-electrons on PMT.
 *
 * The number of photo-electrons of a hit may be larger than unity to limit
 * the overall number of hits and consequently the number of times the arrival
 * time needs to be evaluated which is CPU intensive.
 *
 * \param  npe    number of photo-electrons on PMT
 * \return        number of photo-electrons of hit
 */
inline int getNumberOfPhotoElectrons(const int npe) {
	const int n = npe >> 4;

	if (n == 0)
		return 1;
	if (n >= 32)
		return 32;

	return n;
}

/**
 * Vertex of energy loss of muon.
 */
class JVertex {
public:
	/**
	 * Default constructor.
	 */
	JVertex() :
			z(0.0), t(0.0), E(0.0), Es(0.0) {
	}

	/**
	 * Constructor.
	 *
	 * \param  z       position [m]
	 * \param  t       time     [ns]
	 * \param  E       energy   [GeV]
	 */
	JVertex(const double z, const double t, const double E) {
		this->z = z;
		this->t = t;
		this->E = E;
		this->Es = 0.0;
	}

	/**
	 * Get position.
	 *
	 * \return         position [m]
	 */
	double getZ() const {
		return z;
	}

	/**
	 * Get time.
	 *
	 * \return         time [ns]
	 */
	double getT() const {
		return t;
	}

	/**
	 * Get muon energy.
	 *
	 * \return         energy [GeV]
	 */
	double getE() const {
		return E;
	}

	/**
	 * Get shower energy.
	 *
	 * \return         energy [GeV]
	 */
	double getEs() const {
		return Es;
	}

	/**
	 * Get range of muon.
	 * This method applies only ionisation energy loss.
	 *
	 * \return         range [m]
	 */
	double getRange() const {
		if (E > MASS_MUON * getIndexOfRefraction())
			return (E - MASS_MUON * getIndexOfRefraction())
					/ (gWater.getA() * DENSITY_SEA_WATER);
		else
			return 0.0;
	}

	/**
	 * Get range of muon.
	 *
	 * \param  geane   energy loss
	 * \return         range [m]
	 */
	double getRange(const JGeane& geane) const {
		return geane(E);
	}

	/**
	 * Apply energy loss energy.
	 *
	 * \param  Es      energy [GeV]
	 */
	void applyEloss(const double Es) {
		this->E -= Es;
		this->Es = Es;
	}

	/**
	 * Step.
	 * This method applies only ionisation energy loss.
	 *
	 * \param  ds      step [m]
	 * \return         this vertex
	 */
	JVertex& step(const double ds) {
		z += ds;
		t += ds / getSpeedOfLight();
		E -= ds * gWater.getA() * DENSITY_SEA_WATER;

		return *this;
	}

	/**
	 * Step.
	 *
	 * \param  geane   energy loss
	 * \param  ds      step [m]
	 * \return         this vertex
	 */
	JVertex& step(const JGeane& geane, const double ds) {
		z += ds;
		t += ds / getSpeedOfLight();
		E = geane(E, ds);

		return *this;
	}

protected:
	double z;
	double t;
	double E;
	double Es;
};

/**
 * Muon trajectory.
 */
class JTrack: public std::vector<JVertex> {
public:
	/**
	 * Constructor.
	 *
	 * \param  vertex  muon starting point
	 */
	JTrack(const JVertex& vertex) :
			std::vector<JVertex>(1, vertex) {
	}
};
}

#endif
