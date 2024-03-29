/**
 * \file FiducialVolume.h
 *
 * \ingroup UBXSec
 * 
 * \brief Class def header for a class FiducialVolume
 *
 * @author Marco Del Tutto
 */

/** \addtogroup UBXSec

    @{*/
#ifndef FIDUCIALVOLUME_H
#define FIDUCIALVOLUME_H

#include <iostream>
#include "fhiclcpp/ParameterSet.h"
//#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoAABox.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoVector.h"

namespace ubana{
  
  /**
   \class FiducialVolume
   User defined class FiducialVolume ... these comments are used to generate
   doxygen documentation!
 */

  class FiducialVolume {
    
  public:
    
    /// Default constructor
    FiducialVolume();

    /// Default destructor
    ~FiducialVolume(){}

    /// Configure function parameters
    void Configure(fhicl::ParameterSet const& p, double det_half_height, double det_width, double det_length);

    /// Printd the current configuration
    void PrintConfig();

    /// Returns true if the vertex is in the active volume
    bool VertexInActive(double x, double y, double z);

    /// Returns true if the vertex is in the FV
    bool VertexInFV(double x, double y, double z);
    
    /// Retruns true if the point is contained
    bool PointContain(double x, double y, double z);

    /// Returns true if the vertex is in the FV 
    bool VertexInFV(double* x);

    /// Returns true if the vertex is in the active volume
    bool VertexInActive(TVector3 x);

    /// Returns true if the vertex is in the FV
    bool VertexInFV(TVector3 x);
    
    /// Return true if the point is in the FV
    bool PointContain(TVector3 x);

    /// Returns true if BOTH points are in the FV
    bool TrackContain(TVector3 x1, TVector3 x2);

  protected:

    double _det_half_height;
    double _det_width;
    double _det_length;

    std::vector<double> _border_x_low;
    std::vector<double> _border_x_high;
    std::vector<double> _border_y_low;
    std::vector<double> _border_y_high;
    std::vector<double> _border_z_low;        
    std::vector<double> _border_z_high;

    size_t _n_fv;

    bool _configured = false;
  };
}

#endif
/** @} */ // end of doxygen group 

