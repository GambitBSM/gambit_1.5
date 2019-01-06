#pragma once

#include "gambit/ColliderBit/detectors/BaseDetector.hpp"

namespace Gambit
{

  namespace ColliderBit
  {

    /// A base class for BuckFast simple smearing simulations within ColliderBit.
    template<typename EventT>
    class BuckFast : BaseDetector<EventT>
    {

      public:

        /// Pointers to actual detector response functions
        /// @{
        void(*smearElectronEnergy)(std::vector<HEPUtils::Particle*>&);
        void(*applyElectronEff)(std::vector<HEPUtils::Particle*>&);
        void(*smearMuonMomentum)(std::vector<HEPUtils::Particle*>&);
        void(*applyMuonEff)(std::vector<HEPUtils::Particle*>&);
        void(*smearTaus)(std::vector<HEPUtils::Particle*>&);
        void(*smearJets)(std::vector<HEPUtils::Particle*>&);
        /// @}

        /// @name Event detection simulation.
        ///@{

        bool partonOnly; ///< Chooses between parton only and full event conversion.
        double antiktR; ///< The jet radius used for the anti-kt jet clustering.

        /// A converter for a Pythia8::Event which considers all final state particles.
        /// @note Also performs the jet clustering algorithm.
        void convertParticleEvent(const EventT&, HEPUtils::Event&) const;

        /// A converter for a Pythia8::Event which considers only partonic final states.
        /// @note Also performs the jet clustering algorithm.
        void convertPartonEvent(const EventT&, HEPUtils::Event&) const;

        ///@}

        /// Constructor
        BuckFast() : partonOnly(false)
                   , antiktR(0.4)
                   , smearElectronEnergy(NULL)
                   , applyElectronEff(NULL)
                   , smearMuonMomentum(NULL)
                   , applyMuonEff(NULL)
                   , smearTaus(NULL)
                   , smearJets(NULL)
        {}

        /// Destructor
        virtual ~BuckFast() {}

        /// @name (Re-)Initialization functions
        ///@{

        /// Settings parsing and initialization for any sub-class.
        virtual void init(const std::vector<std::string>&) {};

        /// General init for any collider of this type - no settings version.
        virtual void init() {};

        /// Settings parsing and initialization for sub-classes with parton and jet radius settings only.
        void init(bool parton, double R) { partonOnly=parton; antiktR=R; };

        ///@}

    };

  }
}
