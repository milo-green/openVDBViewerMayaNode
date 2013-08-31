#ifndef __VDBVIEWERNODE__
#define __VDBVIEWERNODE__
#define GL_GLEXT_PROTOTYPES 1

#include <maya/MPxLocatorNode.h>

#include <maya/MVectorArray.h>
#include <maya/MFloatArray.h>
#include <maya/MBoundingBox.h>
#include <maya/MString.h>
#include <maya/MStringArray.h>
#include <maya/MTime.h>
#include <maya/MMatrix.h>
#include <GL/gl.h>
#include <GL/glu.h>


#include <GL/glu.h>
#include <stdlib.h>


#include <openvdb/openvdb.h>
#include <openvdb/util/logging.h>

#include <openvdb/tools/VolumeToMesh.h>

#define USE_CALLBACKS 0
#define USE_CONVERTER 0

typedef std::vector<std::string> StringVec;

class vdbViewerNode : public MPxLocatorNode
{
public:
    static MTypeId id;

    static MObject aFilterChange;
    static MObject aDummy;
    static MObject aVdbFile;
    static MObject aPercent;
    static MObject aChannels;
    static MObject aViewAs;
    static MObject aExposure;
    static MObject aPointSize;
    static MObject aSmooth;
    static MObject aInvert;
    static MObject aAbsoluteValue;
    static MObject aShowVectors;
    static MObject ashowVoxelTree;
    static MObject aVectorsSize;
    static MObject aInvertVectorsColor;
    static MObject aChannelNames;
    static MObject aNumPoints;
    static MObject aNumPointsLoaded;
    static MObject aUseCropBox;
    static MObject aCropBoxMin;
    static MObject aCropBoxMax;
    static MObject aShowProgress;
    static MObject aTime;
    static MObject aShowValue;
    static MObject aFilterVariable;
    static MObject aFilterMode;
    static MObject aFilterValf1;
    static MObject aFilterValf2;
    static MObject aFilterValc1;
    static MObject aFilterValc2;
    static MObject aCircleSlices;
    static MObject aDiskSlices;
    static MObject aFilterOn;

    vdbViewerNode();
    virtual ~vdbViewerNode();

    virtual void postConstructor();

    virtual MStatus compute( const MPlug& plug, MDataBlock& data );

    virtual void draw( M3dView & view,
                    const MDagPath & path,
                    M3dView::DisplayStyle style,
                    M3dView::DisplayStatus displaystatus
                    );

    virtual bool isBounded() const;
    virtual MBoundingBox boundingBox() const;

    static  void* creator();
    static  MStatus initialize();
    MStatus loadVDB();
    double  lerp( double a, double b, double v ) { return (a*(1.0-v)+b*v); };
    void draw_bbox( bool state );
    template<typename GridType>
    void draw_tree(typename GridType::Ptr grid);
    template<typename GridType>
    void draw_mesh(typename GridType::Ptr grid);
    template<typename GridType>
    void draw_meshDummy(typename GridType::Ptr grid);
    void draw_cropBox( bool state );
    void compute_cropBox_world();
    //void compute_cropBox_Vectorized();
    inline bool isInCropBox( float *point );
    void resetDisplayList( int &list );
    void forceUIRefresh();
    MStatus connectDummyToShear();
    MStatus connectToTime();

    inline bool filterPoint( MVector &val, float &radius, bool load );
    void printLongListing();
    void printShortListing( bool metadata);
    void vdbGetBound( bool worldspace, float inbbox[6],MStatus &status);
    void vdbOpenInitilise(const std::string& filename, MStatus &status);
    template<typename GridType>
    void gridValues(typename GridType::Ptr grid);
    //template<typename OpType>
    void processTypedGrid(openvdb::GridBase::Ptr grid);
    void processTypedGridTree(openvdb::GridBase::Ptr grid);
    void processTypedGridMesh(openvdb::GridBase::Ptr grid);
    //void vdbReadDataPoint( MVectorArray &m_pts,MFloatArray &m_ptsRadius,MVectorArray &m_ptsVector,MVectorArray &m_ptsColor,int &m_nPointsLoaded,float s,int numPts, int nPoints, bool doCropBox, bool showProgress, float progressStep);
    void vdbGetNumberNamesTypes( int &numChannels, std::vector<std::string>& names, std::vector<std::string>& types);
    void vdbGetNumberPoints( int &numPoints);
    double fclamp( double x,double a,double b) const {return fmax(a,fmin(x,b));};
    MString buildFileName();



private:
    bool           m_bDummyConnected;
    bool           m_bFilterChangeConnected;
    bool           m_timeConnected;
    MBoundingBox   m_bbox;
    GLUquadricObj *m_quadricObj;
    MString        m_vdbFile;
    MVectorArray   m_pts;
    MVectorArray   m_ptsColor;
    MVectorArray   m_ptsVector;
    MFloatArray    m_ptsRadius;
    MString        m_channelType;
    int            m_chan;
    int            m_viewAs;
    int            m_nPoints;
    int            m_nPointsLoaded;
    float          m_percent;
    float          m_exp;
    float          m_ptSize;
    bool           m_smooth;
    bool           m_invert;
    bool           m_abs;
    bool           m_Vectors;
    bool           m_showVoxelTree;
    float          m_VectorsSize;
    bool           m_invertVectorsColor;
    MStringArray   m_channelNames;
    bool           m_doCropBox;
    MBoundingBox   m_cropBox;
    MBoundingBox   m_cropBoxLocal;
    bool           m_timeSync;
    bool           m_ShowProgress;
    MTime          m_time;
    bool           m_showValue;
    int            m_filterVariable;
    int            m_filterMode;
    float          m_filterValf1;
    float          m_filterValf2;
    MVector        m_filterValc1;
    MVector        m_filterValc2;
    int            m_circleSlices;
    int            m_diskSlices;
    int            m_filterOn;
    // display lists IDs
    int            m_bboxID;
    int            m_treeID;
    int            m_meshID;
    int            m_cropboxID;
    int            m_pcID;
    int            m_diskID;

    std::string m_vdbFileName;
    std::vector<std::string> m_varnames;
    float s;
    int numPts;
    int nPoints;
    bool showProgress;
    float progressStep;
    

    // VDB varibles
    openvdb::GridPtrVecPtr grids;
    openvdb::GridCPtrVec allGrids;
    openvdb::MetaMap::Ptr meta;
    std::string version;

};





#endif // __VDBVIEWERNODE__
