////////////////////////////////////////////////////////////////////////
//
// @file TPCNeutrinoIDAna_module.cc
//
// @brief A basic "skeleton" Ana module to serve as an example/basis
//        for the neutrino ID chain
//
// @authors usher@slac.stanford.edu (cloned from an example)
//
///////////////////////////////////////////////////////////////////////

#ifndef  TPCNeutrinoIDAna_Module
#define  TPCNeutrinoIDAna_Module

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/exception.h"
#include "cetlib/cpu_timer.h"

// LArSoft includes
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/Vertex.h"
#include "AnalysisBase/CosmicTag.h"
#include "Geometry/Geometry.h"
#include "SimpleTypesAndConstants/geo_types.h"
#include "Utilities/TimeService.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <http://root.cern.ch/root/html532/ClassIndex.html>
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TTree.h"

namespace  TPCNeutrinoIDAna
{
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class definition

class  TPCNeutrinoIDAna : public art::EDAnalyzer
{
public:
 
    // Standard constructor and destructor for an ART module.
    explicit  TPCNeutrinoIDAna(fhicl::ParameterSet const& pset);
    virtual ~ TPCNeutrinoIDAna();

    // This method is called once, at the start of the job. In this
    // example, it will define the histograms and n-tuples we'll write.
    void beginJob();

    // This method is called once, at the start of each run. It's a
    // good place to read databases or files that may have
    // run-dependent information.
    void beginRun(const art::Run& run);

    // This method reads in any parameters from the .fcl files. This
    // method is called 'reconfigure' because it might be called in the
    // middle of a job; e.g., if the user changes parameter values in an
    // interactive event display.
    void reconfigure(fhicl::ParameterSet const& pset);

    // The analysis routine, called once per event. 
    void analyze (const art::Event& evt); 

private:

    // Need vectors here because we have have several instantiations
    // of the module running depending on matching
    std::vector<std::string> fVertexModuleLabelVec;
    std::vector<std::string> fVtxTrackAssnsModuleLabelVec;
    
    // Pointers to the histograms we'll create. 

    // The variables that will go into the n-tuple.
    int fEvent;
    int fRun;
    int fSubRun;

    // Other variables that will be shared between different methods.
    art::ServiceHandle<geo::Geometry>            fGeometry;       // pointer to Geometry service
    art::ServiceHandle<util::DetectorProperties> fDetectorProperties;

}; // class  TPCNeutrinoIDAna


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class implementation

//-----------------------------------------------------------------------
// Constructor
 TPCNeutrinoIDAna:: TPCNeutrinoIDAna(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
{
    // Read in the parameters from the .fcl file.
    this->reconfigure(parameterSet);
}

//-----------------------------------------------------------------------
// Destructor
 TPCNeutrinoIDAna::~ TPCNeutrinoIDAna()
{}
   
//-----------------------------------------------------------------------
void  TPCNeutrinoIDAna::beginJob()
{
    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;

    // Define the histograms. Putting semi-colons around the title
    // causes it to be displayed as the x-axis label if the histogram
    // is drawn.
//    fPDGCodeHist        = tfs->make<TH1D>("pdgcodes",";PDG Code;",                  5000, -2500, 2500);
}
   
//-----------------------------------------------------------------------
void  TPCNeutrinoIDAna::beginRun(const art::Run& /*run*/)
{
}

//-----------------------------------------------------------------------
void  TPCNeutrinoIDAna::reconfigure(fhicl::ParameterSet const& pset)
{
    // Read parameters from the .fcl file.
    fVertexModuleLabelVec        = pset.get< std::vector<std::string> >("VertexModuleLabelVec",       std::vector<std::string>() ={"pandoraNu"});
    fVtxTrackAssnsModuleLabelVec = pset.get< std::vector<std::string> >("VtxTrackAssnModuleLabelVec", std::vector<std::string>() ={"neutrinoID"});
    
    if (fVertexModuleLabelVec.size() != fVtxTrackAssnsModuleLabelVec.size())
    {
        mf::LogError("TPCNeutrinoIDFilter") << "Mismatch between string vector lengths input from fhicl!" << std::endl;
    }
    
    return;
}

//-----------------------------------------------------------------------
void  TPCNeutrinoIDAna::analyze(const art::Event& event)
{
    // Start by fetching some basic event information for our n-tuple.
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();
    
    // In principle we can have several producers running over various configurations of vertices and tracks.
    // The output associations we want to check are then encapsuated in the input vectors of strings
    // So the outer loop is over the indices
    for(size_t assnIdx = 0; assnIdx < fVertexModuleLabelVec.size(); assnIdx++)
    {
        // Recover a handle to the collection of vertices
        art::Handle< std::vector<recob::Vertex> > vertexVecHandle;
        event.getByLabel(fVertexModuleLabelVec[assnIdx], vertexVecHandle);
        
        if (vertexVecHandle.isValid())
        {
            // Recover associations relating vertices and tracks
            art::FindManyP<recob::Track> vertexTrackAssns(vertexVecHandle, event, fVtxTrackAssnsModuleLabelVec[assnIdx]);
            
            // First check that we have something
            if (vertexTrackAssns.isValid() && vertexTrackAssns.size() > 0)
            {
                std::cout << ">>>> Vertex/track associations for producer: " << fVtxTrackAssnsModuleLabelVec[assnIdx] << std::endl;
                
                // Loop over vertex/track associations
                for(size_t vtxTrackAssnIdx = 0; vtxTrackAssnIdx < vertexTrackAssns.size(); vtxTrackAssnIdx++)
                {
                    const std::vector<art::Ptr<recob::Track>>& trackVec = vertexTrackAssns.at(vtxTrackAssnIdx);
                    
                    //Loop over tracks
                    for(const auto& track : trackVec)
                    {
                        const TVector3& trackStart = track->Vertex();
                        const TVector3& trackEnd   = track->End();
                        
                        // Geometry routines want to see an array...
                        double trackStartPos[] = {trackStart.X(),trackStart.Y(),trackStart.Z()};
                        double trackEndPos[]   = {trackEnd.X(),  trackEnd.Y(),  trackEnd.Z()};
                        
                        // now loop over views for the construction of the output url
                        for(size_t viewIdx = 0; viewIdx < fGeometry->Nviews(); viewIdx++)
                        {
                            unsigned int startWire  = fGeometry->NearestWire(trackStartPos, viewIdx);
                            double       startTicks = fDetectorProperties->ConvertXToTicks(trackStartPos[0], int(viewIdx), 0, 0);
                            unsigned int endWire    = fGeometry->NearestWire(trackEndPos,   viewIdx);
                            double       endTicks   = fDetectorProperties->ConvertXToTicks(trackEndPos[0], int(viewIdx), 0, 0);
                            
                            // to compile...
                            std::cout << startWire << startTicks << endWire << endTicks << std::endl;
                        }
                    }
                }
            }
        }
    }

    return;
}

DEFINE_ART_MODULE( TPCNeutrinoIDAna)

} // namespace  TPCNeutrinoIDAna

#endif //  TPCNeutrinoIDAna_Module
