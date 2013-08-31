#include "vdbViewerNode.h"

#include <maya/MFnPlugin.h>
#include <maya/MFnStringArrayData.h>
#include <maya/MTypeId.h>
#include <maya/MPlug.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MPoint.h>
#include <maya/MVector.h>
#include <maya/MSelectionList.h>
#include <maya/MDGModifier.h>
#include <maya/MFnDagNode.h>
#include <maya/MDagPath.h>
#include <maya/MAnimControl.h>
#include <maya/MFileObject.h>
#include <maya/MGlobal.h>
#include <maya/M3dView.h>
#include <maya/MProgressWindow.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <openvdb/openvdb.h>


#define vdbViewerNodeId 0x879001


#define VALUE_FILTER  0
#define RADIUS_FILTER 1
#define FILTER_PASSTHROUGH	0
#define FILTER_LESSTHAN		1
#define FILTER_MORETHAN		2
#define FILTER_EQUAL		3
#define FILTER_BETWEEN		4
#define FILTER_ON_LOAD		0
#define FILTER_ON_DRAW		1
#define RAD_TO_DEG 57.29577951308232087679
#define DEG_TO_RAD  .01745329251994329576
#define CHECKERR(STAT,MSG)                                 \
  if ( MS::kSuccess != STAT ) {                            \
	cerr <<"  [ Failed ] : " <<MSG << " : "<<STAT<< endl;  \
	return MS::kFailure;                                   \
  }


#ifdef _DEBUG
#   define DEBUG(msg)    cout <<"   [ vdbViewerNode ]  "<<msg<<endl<<flush
#   define DEBUGV(msg,x) cout <<"   [ vdbViewerNode ]  "<<msg<<x<<endl<<flush
#else
#   define DEBUG(msg)
#   define DEBUGV(msg,x)
#endif

#define ERR(msg)                                        \
cerr <<"   [ vdbViewerNode ]  ERROR : "<<msg<<endl;   \
MGlobal::displayError( msg );

#define ERRV(msg,x)                                        \
cerr <<"   [ vdbViewerNode ]  ERROR : "<<msg<<x<<endl;   \
MGlobal::displayError( msg );


MTypeId vdbViewerNode::id( vdbViewerNodeId );

MObject vdbViewerNode::aFilterChange;
MObject vdbViewerNode::aDummy;
MObject vdbViewerNode::aVdbFile;
MObject vdbViewerNode::aChannels;
MObject vdbViewerNode::aViewAs;
MObject vdbViewerNode::aPercent;
MObject vdbViewerNode::aExposure;
MObject vdbViewerNode::aPointSize;
MObject vdbViewerNode::aSmooth;
MObject vdbViewerNode::aInvert;
MObject vdbViewerNode::aAbsoluteValue;
MObject vdbViewerNode::aShowVectors;
MObject vdbViewerNode::ashowVoxelTree;
MObject vdbViewerNode::aVectorsSize;
MObject vdbViewerNode::aInvertVectorsColor;
MObject vdbViewerNode::aChannelNames;
MObject vdbViewerNode::aNumPoints;
MObject vdbViewerNode::aNumPointsLoaded;
MObject vdbViewerNode::aUseCropBox;
MObject vdbViewerNode::aCropBoxMin;
MObject vdbViewerNode::aCropBoxMax;
MObject vdbViewerNode::aShowProgress;
MObject vdbViewerNode::aTime;
MObject vdbViewerNode::aShowValue;
MObject vdbViewerNode::aFilterVariable;
MObject vdbViewerNode::aFilterMode;
MObject vdbViewerNode::aFilterValf1;
MObject vdbViewerNode::aFilterValf2;
MObject vdbViewerNode::aFilterValc1;
MObject vdbViewerNode::aFilterValc2;
MObject vdbViewerNode::aCircleSlices;
MObject vdbViewerNode::aDiskSlices;
MObject vdbViewerNode::aFilterOn;

typedef std::vector<std::string> StringVec;

const char* INDENT = "   ";
const char* INDENT2 = "       ";

std::string
sizeAsString(openvdb::Index64 n, const std::string& units)
{
    std::ostringstream ostr;
    ostr << std::setprecision(3);
    if (n < 1000) {
        ostr << n;
    } else if (n < 1000000) {
        ostr << (n / 1.0e3) << "K";
    } else if (n < 1000000000) {
        ostr << (n / 1.0e6) << "M";
    } else {
        ostr << (n / 1.0e9) << "G";
    }
    ostr << units;
    return ostr.str();
}


std::string
bytesAsString(openvdb::Index64 n)
{
    std::ostringstream ostr;
    ostr << std::setprecision(3);
    if (n >> 30) {
        ostr << (n / double(uint64_t(1) << 30)) << "GB";
    } else if (n >> 20) {
        ostr << (n / double(uint64_t(1) << 20)) << "MB";
    } else if (n >> 10) {
        ostr << (n / double(uint64_t(1) << 10)) << "KB";
    } else {
        ostr << n << "B";
    }
    return ostr.str();
}


std::string
coordAsString(const openvdb::Coord ijk, const std::string& sep)
{
    std::ostringstream ostr;
    ostr << ijk[0] << sep << ijk[1] << sep << ijk[2];
    return ostr.str();
}


/// Return a string representation of the given metadata key, value pairs
std::string
metadataAsString(
    const openvdb::MetaMap::ConstMetaIterator& begin,
    const openvdb::MetaMap::ConstMetaIterator& end,
    const std::string& indent = "")
{
    std::ostringstream ostr;
    char sep[2] = { 0, 0 };
    for (openvdb::MetaMap::ConstMetaIterator it = begin; it != end; ++it) {
        ostr << sep << indent << it->first;
        if (it->second) {
            const std::string value = it->second->str();
            if (!value.empty()) ostr << ": " << value;
        }
        sep[0] = '\n';
    }
    return ostr.str();
}


std::string
bkgdValueAsString(const openvdb::GridBase::ConstPtr& grid)
{
    std::ostringstream ostr;
    if (grid) {
        const openvdb::TreeBase& tree = grid->baseTree();
        ostr << "background: ";
        openvdb::Metadata::Ptr background = tree.getBackgroundValue();
        if (background) ostr << background->str();
    }
    return ostr.str();
}


/// Print detailed information about the given VDB files.
/// If @a metadata is true, include file-level metadata key, value pairs.
void vdbViewerNode::printLongListing()
{

        // Print file-level metadata.
        std::cout << "VDB version: " << version << "\n";
        if (meta) {
            std::string str = metadataAsString(meta->beginMeta(), meta->endMeta());
            if (!str.empty()) std::cout << str << "\n";
        }
        std::cout << "\n";

        // For each grid in the file...
        bool firstGrid = true;
        for (openvdb::GridPtrVec::const_iterator it = grids->begin(); it != grids->end(); ++it) {
            if (openvdb::GridBase::ConstPtr grid = *it) {
                if (!firstGrid) std::cout << "\n\n";
                std::cout << "Name: " << grid->getName() << std::endl;
                grid->print(std::cout, /*verboseLevel=*/3);
                firstGrid = false;
            }
        }

}


/// Print condensed information about the given VDB files.
/// If @a metadata is true, include file- and grid-level metadata.
void vdbViewerNode::printShortListing( bool metadata)
{
        std::cout <<"\n------------------------------\n";
        std::cout << INDENT << std::left << std::setw(11) << "Vdb File: "<< m_vdbFileName << "\n";
        std::cout <<"------------------------------\n";

        if (!grids);

        if (metadata) {
            // Print file-level metadata.
            std::string str = metadataAsString(meta->beginMeta(), meta->endMeta(), INDENT);
            if (!str.empty()) std::cout << str << "\n";
        }

        // For each grid in the file...
        for (openvdb::GridPtrVec::const_iterator it = grids->begin(); it != grids->end(); ++it) {
            const openvdb::GridBase::ConstPtr grid = *it;
            if (!grid) continue;

            // Print the grid name and its voxel value datatype.
            std::cout << INDENT << std::left << std::setw(11) << grid->getName()
                << " " << std::right << std::setw(6) << grid->valueType();

            // Print the grid's bounding box and dimensions.
            openvdb::CoordBBox bbox = grid->evalActiveVoxelBoundingBox();
            std::string
                boxStr = coordAsString(bbox.min()," ") + "  " + coordAsString(bbox.max()," "),
                dimStr = coordAsString(bbox.extents(), "x");
            boxStr += std::string(std::max<int>(1,
                40 - boxStr.size() - dimStr.size()), ' ') + dimStr;
            std::cout << " " << std::left << std::setw(40) << boxStr;

            // Print the number of active voxels.
            std::cout << "  " << std::right << std::setw(8)
                << sizeAsString(grid->activeVoxelCount(), "Vox");

            // Print the grid's in-core size, in bytes.
            std::cout << " " << std::right << std::setw(6) << bytesAsString(grid->memUsage());

            std::cout << std::endl;

            // Print grid-specific metadata.
            if (metadata) {
                // Print background value.
                std::string str = bkgdValueAsString(grid);
                if (!str.empty()) {
                    std::cout << INDENT2 << str << "\n";
                }
                // Print local and world transforms.
                grid->transform().print(std::cout, INDENT2);
                // Print custom metadata.
                str = metadataAsString(grid->beginMeta(), grid->endMeta(), INDENT2);
                if (!str.empty()) std::cout << str << "\n";
                std::cout << std::flush;
            }
        }

}


/// get the bound
void vdbViewerNode::vdbGetBound( bool worldspace, float inbbox[6],MStatus &status)
{

        if (grids->empty()) {
                //OPENVDB_LOG_WARN( " grids are empty");
                status =  MS::kFailure;
                return;
            }
        openvdb::GridCPtrVec allGrids;
        allGrids.insert(allGrids.end(), grids->begin(), grids->end());

        // itterate grids find min max bund box
        openvdb::Vec3d min(std::numeric_limits<double>::max()), max(-min);
        for (size_t n = 0; n < allGrids.size(); ++n) {
        openvdb::CoordBBox bbox = allGrids[n]->evalActiveVoxelBoundingBox();
        min = openvdb::math::minComponent(min, allGrids[n]->indexToWorld(bbox.min()));
        max = openvdb::math::maxComponent(max, allGrids[n]->indexToWorld(bbox.max()));
        }

        inbbox[0] = min.x();
        inbbox[1] = min.y();
        inbbox[2] = min.z();
        inbbox[3] = max.x();
        inbbox[4] = max.y();
        inbbox[5] = max.z();
}

void vdbViewerNode::vdbOpenInitilise(const std::string& filename, MStatus &status)
{
openvdb::io::File file(filename);
        try {
            file.open();
            grids = file.getGrids();
            meta = file.getMetadata();
            version = file.version();
            allGrids.clear();
            allGrids.insert(allGrids.end(), grids->begin(), grids->end());
            file.close();
        } catch (openvdb::Exception& e) {
            OPENVDB_LOG_ERROR(e.what() << " (" << filename << ")");
            status =  MS::kFailure;
            return;
        }
        if (grids->empty()) {
                OPENVDB_LOG_WARN(filename << " is empty");
                return;
            }
        return;
}





void vdbViewerNode::vdbGetNumberNamesTypes( int &numChannels, std::vector<std::string>& names, std::vector<std::string>& types)
{

        // For each grid in the file... list name and type
        for (openvdb::GridPtrVec::const_iterator it = grids->begin(); it != grids->end(); ++it) {
            const openvdb::GridBase::ConstPtr grid = *it;
            if (!grid) continue;

            numChannels+=1;
            names.push_back(grid->getName());
            types.push_back(grid->valueType());

            
        }
    
}

void vdbViewerNode::vdbGetNumberPoints( int &numPoints)
{
        // For each grid in the file...
        for (openvdb::GridPtrVec::const_iterator it = grids->begin(); it != grids->end(); ++it) {
            const openvdb::GridBase::ConstPtr grid = *it;
            if (!grid) continue;

            numPoints= grid->activeVoxelCount();
            
        }
}




template< typename valueType>
void processValue(valueType val,MVector &Colour )
{
     Colour = MVector( double(val),double(val),double(val));
}


template<>
void processValue( openvdb::Vec3i val,MVector &Colour )
{
    Colour  = MVector( double(val.x()),double(val.y()),double(val.z()));
}

template<>
void processValue( openvdb::Vec3d val,MVector &Colour )
{
    Colour  = MVector( double(val.x()),double(val.y()),double(val.z()));
}


template<>
void processValue( openvdb::Vec3f val,MVector &Colour )
{
    Colour  = MVector( double(val.x()),double(val.y()),double(val.z()));
}


template<typename GridType>
    void processScalarType(GridType val,MVector &Colour ) 
    {
        //typename GridType  value = it.getValue();
        Colour = MVector( double(val),double(val),double(val));

    }

template<typename GridType>
    void processVectorType(GridType val,MVector &Colour ) 
    {
        //typename GridType value = it.getValue();
        Colour  = MVector( double(val.x()),double(val.y()),double(val.z()));

    }

template<typename GridType>
    void vdbViewerNode::gridValues(typename GridType::Ptr theGrid) 
    {
        int i =0;
        int sp = 0;
        float point[3];
        MVector faraway = MVector(1000000.0,1000000.0,1000000.0);
        float radius;

        for (typename GridType::ValueOnCIter iter = theGrid->cbeginValueOn(); iter; ++iter) 
        {
            openvdb::Vec3d worldSpace  = theGrid->indexToWorld(iter.getCoord());


                            m_pts.set( faraway, sp );

                            point[0] = worldSpace.x();
                            point[1] = worldSpace.y();
                            point[2] = worldSpace.z();
                            
                            bool insideROI = true;
                            //if inCropBox is on test and skip if nessasarry
                            if(m_doCropBox)
                            {
                                if ( !isInCropBox( point ) )
                                {
                                    insideROI=false;
                                }
                            }

                            if(insideROI)
                            {
                                if ( fmod((float)i,s) == 0.0 )
                                {
                                    
                                    if ( sp < numPts )
                                    {


                                        MVector   pColor;
                                        processValue(iter.getValue(), pColor );


                                        radius = 0.05 + pColor.length()*0.01;
                                        if ( filterPoint(pColor,radius,true) )
                                        {

                                            
                                            m_pts.set( point, sp );                                            
                                            m_ptsRadius.set( radius, sp );
                                            m_ptsVector.set( pColor, sp );
                                            m_ptsColor.set( pColor, sp );
                                            m_nPointsLoaded++;
                                        }
                                    }
                                    sp++;
                                }
                            }   
                    i++;
                    if ( showProgress && ( fmod((float)i,progressStep) == 0.0 || i >= nPoints ) )  MProgressWindow:: setProgress(i);
        }
    }



void vdbViewerNode::processTypedGrid(openvdb::GridBase::Ptr grid)
{
//#define CALL_OP(GridType) \
    //op.template operator()<GridType>(openvdb::gridPtrCast<GridType>(grid))

#define CALL_OP(GridType) \
    gridValues<GridType>(openvdb::gridPtrCast<GridType>(grid))

    if (grid->isType<openvdb::BoolGrid>())        CALL_OP(openvdb::BoolGrid);
    else if (grid->isType<openvdb::FloatGrid>())  CALL_OP(openvdb::FloatGrid);
    else if (grid->isType<openvdb::DoubleGrid>()) CALL_OP(openvdb::DoubleGrid);
    else if (grid->isType<openvdb::Int32Grid>())  CALL_OP(openvdb::Int32Grid);
    else if (grid->isType<openvdb::Int64Grid>())  CALL_OP(openvdb::Int64Grid);
    else if (grid->isType<openvdb::Vec3IGrid>())  CALL_OP(openvdb::Vec3IGrid);
    else if (grid->isType<openvdb::Vec3SGrid>())  CALL_OP(openvdb::Vec3SGrid);
    else if (grid->isType<openvdb::Vec3DGrid>())  CALL_OP(openvdb::Vec3DGrid);
    //need to return failure here if type unknow i.e string grid

#undef CALL_OP
}

void vdbViewerNode::processTypedGridTree(openvdb::GridBase::Ptr grid)
{
//#define CALL_OP(GridType) \
    //op.template operator()<GridType>(openvdb::gridPtrCast<GridType>(grid))

#define CALL_OP(GridType) \
    draw_tree<GridType>(openvdb::gridPtrCast<GridType>(grid))

    if (grid->isType<openvdb::BoolGrid>())        CALL_OP(openvdb::BoolGrid);
    else if (grid->isType<openvdb::FloatGrid>())  CALL_OP(openvdb::FloatGrid);
    else if (grid->isType<openvdb::DoubleGrid>()) CALL_OP(openvdb::DoubleGrid);
    else if (grid->isType<openvdb::Int32Grid>())  CALL_OP(openvdb::Int32Grid);
    else if (grid->isType<openvdb::Int64Grid>())  CALL_OP(openvdb::Int64Grid);
    else if (grid->isType<openvdb::Vec3IGrid>())  CALL_OP(openvdb::Vec3IGrid);
    else if (grid->isType<openvdb::Vec3SGrid>())  CALL_OP(openvdb::Vec3SGrid);
    else if (grid->isType<openvdb::Vec3DGrid>())  CALL_OP(openvdb::Vec3DGrid);
    //need to return failure here if type unknow i.e string grid

#undef CALL_OP
}


void vdbViewerNode::processTypedGridMesh(openvdb::GridBase::Ptr grid)
{
//#define CALL_OP(GridType) \
    //op.template operator()<GridType>(openvdb::gridPtrCast<GridType>(grid))

#define CALL_OP(GridType) \
    draw_mesh<GridType>(openvdb::gridPtrCast<GridType>(grid))

    #define CALL_OPDUMMY(GridType) \
    draw_meshDummy<GridType>(openvdb::gridPtrCast<GridType>(grid))

    if (grid->isType<openvdb::BoolGrid>())        CALL_OPDUMMY(openvdb::BoolGrid);
    else if (grid->isType<openvdb::FloatGrid>())  CALL_OP(openvdb::FloatGrid);
    else if (grid->isType<openvdb::DoubleGrid>()) CALL_OP(openvdb::DoubleGrid);
    else if (grid->isType<openvdb::Int32Grid>())  CALL_OPDUMMY(openvdb::Int32Grid);
    else if (grid->isType<openvdb::Int64Grid>())  CALL_OPDUMMY(openvdb::Int64Grid);
    else if (grid->isType<openvdb::Vec3IGrid>())  CALL_OPDUMMY(openvdb::Vec3IGrid);
    else if (grid->isType<openvdb::Vec3SGrid>())  CALL_OPDUMMY(openvdb::Vec3SGrid);
    else if (grid->isType<openvdb::Vec3DGrid>())  CALL_OPDUMMY(openvdb::Vec3DGrid);
    //need to return failure here if type unknow i.e string grid

#undef CALL_OP
}




MString vdbViewerNode::buildFileName()
{
    //string to hold filename
    
    MString fileName = m_vdbFile;
    const float frame = m_time.value();
    
    MString newFileName;
    //int doPadding = 0;
    int padding = 0;
    
    //check if it contains a hash
    int hasDollarF = fileName.indexW(MString("$F"));
    int isFloatFrame  = fileName.indexW(MString("$FF"));

    MString sFrame = MString("");
    if(isFloatFrame!=-1)
    {
        float f;
        f = float ( floor(frame*100+0.5)/100.0 );
        sFrame.set(f,2);
    }
    else
    {
        sFrame+=int(frame);
    }

    int dollarFindex = 0;
    if(hasDollarF != -1)
    {
        m_timeSync = true;
        MStringArray tempStringArray;
        fileName.split ('.',  tempStringArray);

        for(unsigned int split =0; split <=tempStringArray.length()-1;split++)
        {
            int found  = tempStringArray[split].indexW(MString("$F"));
            if(found != -1)
            {
                dollarFindex =split;
            }

        }
        MString frameString = tempStringArray[dollarFindex];

        //if the length is 3  and cant find $FF we assume they have put $F4 or $F3 for padding
        if((frameString.length()==3) && (isFloatFrame==-1))
        {
            MString paddingChar  = frameString.substring(2,3);
            padding = paddingChar.asInt();
            int numZeros = padding - sFrame.length();
            for(int j =0; j<numZeros;j++)
            {
                sFrame = MString("0") + sFrame;
            }
            
        }


        tempStringArray[dollarFindex] = sFrame;
        newFileName = "";
        for(unsigned int split =0; split <=tempStringArray.length()-1;split++)
        {
            if(split!=0)
            {
                newFileName += MString(".");
            }
            newFileName  += tempStringArray[split];
            

        }
        
    
    }
    else
    {
    newFileName = fileName;
    }
    //MGlobal::displayInfo( newFileName );
    
    return newFileName ;
}


//
//  vdbViewerNode::Constructor
//
vdbViewerNode::vdbViewerNode()
{
	m_pts.clear();
	m_ptsColor.clear();
	m_ptsVector.clear();
	m_ptsRadius.clear();
	m_channelNames.clear();
	m_bDummyConnected     = false;
	m_bFilterChangeConnected = false;
	m_timeConnected       = false;
	m_quadricObj          = gluNewQuadric();
	m_vdbFile             = "";
	m_nPoints             = 0;
	m_nPointsLoaded       = 0;
	m_percent             = 10.0;
	m_chan                = 0;
	m_viewAs              = 0;
	m_channelType         = "";
	m_exp                 = 0.0;
	m_ptSize              = 0.0;
	m_smooth              = false;
	m_invert              = false;
	m_abs                 = false;
	m_Vectors             = false;
	m_VectorsSize         = 1.0;
	m_invertVectorsColor  = false;
	m_bbox                = MBoundingBox( MPoint( -1.0, -1.0, -1.0 ), MPoint( 1.0, 1.0, 1.0 ) );
	m_doCropBox           = false;
	m_cropBox             = MBoundingBox( MPoint(0.5,0.0,0.0), MPoint(0.55,1.0,1.0) );
	m_cropBoxLocal        = MBoundingBox( MPoint(0.5,0.0,0.0), MPoint(0.55,1.0,1.0) );
	m_timeSync            = false;
	m_time.setValue(0.0);
	m_showValue           = false;
	m_circleSlices        = 8;
	m_diskSlices          = 8;
	m_filterOn            = 0;
	m_filterVariable      = 1;
	m_filterMode          = 0;
	m_filterValf1         = 0.0;
	m_filterValf2         = 0.0;
	m_filterValc1         = MVector(0.0,0.0,0.0);
	m_filterValc2         = MVector(1.0,1.0,1.0);;

	m_bboxID    = -1;
	m_cropboxID = -1;
	m_pcID     = -1;
	m_diskID    = -1;
    m_treeID    = -1;
    m_meshID    = -1;
}


//
//  vdbViewerNode::Destructor
//
vdbViewerNode::~vdbViewerNode()
{
	resetDisplayList( m_bboxID );
    resetDisplayList( m_treeID );
    resetDisplayList( m_meshID );
	resetDisplayList( m_cropboxID );
	resetDisplayList( m_pcID );
	resetDisplayList( m_diskID );

	if(m_quadricObj == NULL) gluDeleteQuadric(m_quadricObj);
}


void vdbViewerNode::postConstructor()
{
  MFnDependencyNode nodeFn(thisMObject());
  nodeFn.setName("vdbViewerShape#");
}


//
//  vdbViewerNode::compute_cropBox_world
//
void vdbViewerNode::compute_cropBox_world()
{
	m_cropBoxLocal.clear();
	m_cropBoxLocal = MBoundingBox(
						MPoint( fclamp( lerp( m_bbox.min().x, m_bbox.max().x, (double)m_cropBox.min().x ), m_bbox.min().x, m_bbox.max().x),
								fclamp( lerp( m_bbox.min().y, m_bbox.max().y, (double)m_cropBox.min().y ), m_bbox.min().y, m_bbox.max().y),
								fclamp( lerp( m_bbox.min().z, m_bbox.max().z, (double)m_cropBox.min().z ), m_bbox.min().z, m_bbox.max().z) ),
						MPoint( fclamp( lerp( m_bbox.min().x, m_bbox.max().x, (double)m_cropBox.max().x ), m_bbox.min().x, m_bbox.max().x),
								fclamp( lerp( m_bbox.min().y, m_bbox.max().y, (double)m_cropBox.max().y ), m_bbox.min().y, m_bbox.max().y),
								fclamp( lerp( m_bbox.min().z, m_bbox.max().z, (double)m_cropBox.max().z ), m_bbox.min().z, m_bbox.max().z) )
								);

}



//
//  vdbViewerNode::isInCropBox
//
inline bool vdbViewerNode::isInCropBox( float *point )
{
	if (    point[0] > m_cropBoxLocal.min().x && point[0] < m_cropBoxLocal.max().x &&
			point[1] > m_cropBoxLocal.min().y && point[1] < m_cropBoxLocal.max().y &&
			point[2] > m_cropBoxLocal.min().z && point[2] < m_cropBoxLocal.max().z
	) return true;
	return false;
}


//
//  vdbViewerNode::resetDisplayList
//
void vdbViewerNode::resetDisplayList( int &list )
{
	if ( list == -1 ) return;
	glDeleteLists( list, 1 );
	list = -1;
	DEBUG(" - Reset list");
}


//
//  vdbViewerNode::forceUIRefresh
//
void vdbViewerNode::forceUIRefresh()
{
	if ( MGlobal::mayaState() == MGlobal::kInteractive ) {
		DEBUG( "REFRESH AE" );
		MGlobal::executeCommandOnIdle( "vdbForceAEUpdate();", false );
	}
}


//
//  vdbViewerNode::compute
//
MStatus vdbViewerNode::compute( const MPlug& plug, MDataBlock& data )
{
	bool reload   = false;
	bool refilter = false;
	bool redraw   = false;
	MStatus returnStatus = MS::kUnknownParameter;
	MStatus status;

	DEBUG("COMPUTE --------------------------------------------");

	if ( plug == aFilterChange )
	{
		DEBUG(" + filter change");

		int old_filterOn = m_filterOn;

		m_filterVariable      = data.inputValue( aFilterVariable ).asShort();
		m_filterMode          = data.inputValue( aFilterMode ).asShort();
		m_filterValf1         = data.inputValue( aFilterValf1 ).asFloat();
		m_filterValf2         = data.inputValue( aFilterValf2 ).asFloat();
		m_filterValc1         = data.inputValue( aFilterValc1 ).asFloatVector();
		m_filterValc2         = data.inputValue( aFilterValc2 ).asFloatVector();
		m_filterOn            = data.inputValue( aFilterOn ).asShort();

		DEBUGV("   + m_filterOn = ",m_filterOn);
		DEBUGV("   + filter color 1: ",m_filterValc1);

		// update to do list
		//
		refilter = true;
		if ( m_filterOn == FILTER_ON_LOAD ) reload = true;
		if ( old_filterOn != m_filterOn ) reload = true;

		returnStatus = MS::kSuccess;
	}

	if ( plug == aDummy )
	{
		DEBUG(" + param change");

		MString old_vdbFile       = m_vdbFile;
		int     old_channel       = m_chan;
		bool    old_doCropBox     = m_doCropBox;
		float   old_percent       = m_percent;
		MTime   old_time          = m_time;
		bool    old_timeSync      = m_timeSync;
		MBoundingBox old_cropBox  = m_cropBox;

		// get input values
		//
		m_vdbFile             = data.inputValue( aVdbFile ).asString();
		m_exp                 = powf( 2.0f, data.inputValue( aExposure ).asFloat() );
		m_ptSize              = data.inputValue( aPointSize ).asFloat();
		m_smooth              = data.inputValue( aSmooth ).asBool();
		m_percent             = data.inputValue( aPercent ).asFloat();
		m_chan                = data.inputValue( aChannels ).asInt();
		m_viewAs              = data.inputValue( aViewAs ).asInt();
		m_invert              = data.inputValue( aInvert ).asBool();
		m_abs                 = data.inputValue( aAbsoluteValue ).asBool();
		m_Vectors             = data.inputValue( aShowVectors ).asBool();
        m_showVoxelTree       = data.inputValue( ashowVoxelTree ).asBool();
		m_VectorsSize         = data.inputValue( aVectorsSize ).asFloat();
		m_invertVectorsColor  = data.inputValue( aInvertVectorsColor ).asBool();

		// Cosmetic adjustment
		// Is smoot points are on, increase the point size by one
		// to keep the perceived point size constant.
		//
		if ( m_smooth ) m_ptSize += 1;

		// get the crop box values
		// we will limit the values to make sure it
		// remains inside the main bbox and that the crop box
		// doesn't fold over itself.
		//
		MPoint minP           = data.inputValue( aCropBoxMin ).asFloatVector();
		MPoint maxP           = data.inputValue( aCropBoxMax ).asFloatVector();
		minP.x = fmin(minP.x, m_cropBox.max().x);
		minP.y = fmin(minP.y, m_cropBox.max().y);
		minP.z = fmin(minP.z, m_cropBox.max().z);
		maxP.x = fmax(maxP.x, m_cropBox.min().x);
		maxP.y = fmax(maxP.y, m_cropBox.min().y);
		maxP.z = fmax(maxP.z, m_cropBox.min().z);
		// make sure the attributes are appropriately limited
		if ( minP != m_cropBox.min() )
		{
			MFnNumericAttribute minAttrFn(aCropBoxMin);
			minAttrFn.setMax( fmin((double)maxP.x, m_bbox.max().x),
							  fmin((double)maxP.y, m_bbox.max().y),
							  fmin((double)maxP.z, m_bbox.max().z) );
		} else if ( maxP != m_cropBox.max() )
		{
			MFnNumericAttribute maxAttrFn(aCropBoxMax);
			maxAttrFn.setMin( fmin((double)minP.x, m_bbox.min().x),
							  fmin((double)minP.y, m_bbox.min().y),
							  fmin((double)minP.z, m_bbox.min().z) );
		}
		m_cropBox             = MBoundingBox(minP,maxP);
		m_doCropBox           = data.inputValue( aUseCropBox ).asBool();

		// the rest of the inputs
		//
		m_ShowProgress            = data.inputValue( aShowProgress ).asBool();
		m_time                = data.inputValue( aTime ).asTime();
		m_showValue           = data.inputValue( aShowValue ).asBool();
		m_circleSlices        = data.inputValue( aCircleSlices ).asInt();
		m_diskSlices          = data.inputValue( aDiskSlices ).asInt();


		//
		// update to do list
		//
		redraw = true;

		if ( reload == true ||
			 old_vdbFile   !=  m_vdbFile  ||
			 old_channel   !=  m_chan     ||
			 old_percent   !=  m_percent  ||
			 old_timeSync  !=  m_timeSync ||
			 ( m_timeSync && old_time != m_time ) ||
			 old_doCropBox !=  m_doCropBox ||
			 ( m_doCropBox && ( old_cropBox.min() != m_cropBox.min() || old_cropBox.max() != m_cropBox.max() ) )
		   ) reload = true;

		returnStatus = MS::kSuccess;
	}

	DEBUGV("   + reload = ",reload);

	if ( reload )
	{
		if ( m_vdbFile != "" )
		{

			MFileObject fileObj;
			fileObj.setRawFullName(buildFileName());

			if ( ! fileObj.exists() )
			{

				// file does not exist
				m_vdbFile = "";
				m_pts.clear();
				m_ptsColor.clear();
				m_ptsVector.clear();
				m_ptsRadius.clear();
				ERR("VDB file does not exist..!");

			}
			else
			{
				DEBUG(" + reload vdb");

				status = loadVDB();
				resetDisplayList( m_bboxID );
                resetDisplayList( m_treeID );
                resetDisplayList( m_meshID );
                resetDisplayList( m_pcID );
                resetDisplayList( m_cropboxID );
                resetDisplayList( m_diskID );
				forceUIRefresh();

				if ( status != MS::kSuccess )
				{
					return MS::kFailure;
				}

				//
				// fill output attributes
				//

				// set the name of the channels
				//
				MDataHandle namesHdl = data.outputValue( aChannelNames, &status );
				if ( status != MS::kSuccess )
				{
					ERR("getting aChannelNames array handle");
					return MS::kFailure;
				}
				MFnStringArrayData strData;
				MObject strObj = strData.create( m_channelNames );
				status = namesHdl.set( strObj );
				if ( status != MS::kSuccess )
				{
					ERR("writing channelNames to attr");
					return MS::kFailure;
				}
				namesHdl.setClean();

				// set the number of points
				//
				MDataHandle numHdl = data.outputValue( aNumPoints, &status );
				numHdl.set( m_nPoints );
				numHdl.setClean();

				// set the number of points loaded
				//
				numHdl = data.outputValue( aNumPointsLoaded, &status );
				numHdl.set( m_nPointsLoaded );
				numHdl.setClean();
			}
		}
	}

	// safer to recompute every time.
	//
	compute_cropBox_world();
	resetDisplayList( m_bboxID );
    resetDisplayList( m_treeID );
    resetDisplayList( m_meshID );
	resetDisplayList( m_cropboxID );
	resetDisplayList( m_pcID );
	resetDisplayList( m_diskID );
	return returnStatus;
}


//
//  vdbViewerNode::draw_bbox
//
void vdbViewerNode::draw_bbox( bool state )
{
	if ( m_bboxID < 0 )
	{

		m_bboxID = glGenLists(1);
		glNewList( m_bboxID, GL_COMPILE_AND_EXECUTE );

			glPushAttrib( GL_ALL_ATTRIB_BITS);

			if ( state == true )
			{
				short pattern = 37449;
				glLineStipple( 1, pattern );
				glEnable( GL_LINE_STIPPLE );
			}

			if ( m_smooth )
			{
				glEnable(GL_BLEND);
				glEnable(GL_LINE_SMOOTH);
				glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
			}

			glBegin( GL_LINES );

				glVertex3f( m_bbox.min().x, m_bbox.max().y, m_bbox.max().z );
				glVertex3f( m_bbox.max().x, m_bbox.max().y, m_bbox.max().z );

				glVertex3f( m_bbox.min().x, m_bbox.max().y, m_bbox.min().z );
				glVertex3f( m_bbox.max().x, m_bbox.max().y, m_bbox.min().z );

				glVertex3f( m_bbox.min().x, m_bbox.min().y, m_bbox.max().z );
				glVertex3f( m_bbox.max().x, m_bbox.min().y, m_bbox.max().z );

				glVertex3f( m_bbox.min().x, m_bbox.min().y, m_bbox.min().z );
				glVertex3f( m_bbox.max().x, m_bbox.min().y, m_bbox.min().z );

				glVertex3f( m_bbox.min().x, m_bbox.min().y, m_bbox.min().z );
				glVertex3f( m_bbox.min().x, m_bbox.min().y, m_bbox.max().z );

				glVertex3f( m_bbox.min().x, m_bbox.max().y, m_bbox.min().z );
				glVertex3f( m_bbox.min().x, m_bbox.max().y, m_bbox.max().z );

				glVertex3f( m_bbox.max().x, m_bbox.min().y, m_bbox.min().z );
				glVertex3f( m_bbox.max().x, m_bbox.min().y, m_bbox.max().z );

				glVertex3f( m_bbox.max().x, m_bbox.max().y, m_bbox.min().z );
				glVertex3f( m_bbox.max().x, m_bbox.max().y, m_bbox.max().z );

				glVertex3f( m_bbox.min().x, m_bbox.min().y, m_bbox.min().z );
				glVertex3f( m_bbox.min().x, m_bbox.max().y, m_bbox.min().z );

				glVertex3f( m_bbox.max().x, m_bbox.min().y, m_bbox.min().z );
				glVertex3f( m_bbox.max().x, m_bbox.max().y, m_bbox.min().z );

				glVertex3f( m_bbox.max().x, m_bbox.min().y, m_bbox.max().z );
				glVertex3f( m_bbox.max().x, m_bbox.max().y, m_bbox.max().z );

				glVertex3f( m_bbox.min().x, m_bbox.min().y, m_bbox.max().z );
				glVertex3f( m_bbox.min().x, m_bbox.max().y, m_bbox.max().z );

			glEnd();

			// center coord sys
			float length = fmin(fmin( m_bbox.width(), m_bbox.height() ), m_bbox.depth() ) * 0.1;
			glDisable( GL_LINE_STIPPLE );
			glBegin( GL_LINES );

				glColor3f( 1.0, 0.0, 0.0 );
				glVertex3f( m_bbox.center().x, m_bbox.center().y, m_bbox.center().z );
				glVertex3f( m_bbox.center().x+length, m_bbox.center().y, m_bbox.center().z );

				glColor3f( 0.0, 1.0, 0.0 );
				glVertex3f( m_bbox.center().x, m_bbox.center().y, m_bbox.center().z );
				glVertex3f( m_bbox.center().x, m_bbox.center().y+length, m_bbox.center().z );

				glColor3f( 0.0, 0.0, 1.0 );
				glVertex3f( m_bbox.center().x, m_bbox.center().y, m_bbox.center().z );
				glVertex3f( m_bbox.center().x, m_bbox.center().y, m_bbox.center().z+length );

			glEnd();

			glPopAttrib();

		glEndList();

	}
	else
	{
		glCallList( m_bboxID );
	}

}

template<typename GridType>
void vdbViewerNode::draw_tree(typename GridType::Ptr grid)
{
    if ( m_treeID < 0 )
    {
        bool state = true;

        m_treeID = glGenLists(1);
        glNewList( m_treeID, GL_COMPILE_AND_EXECUTE );

            glPushAttrib( GL_ALL_ATTRIB_BITS);

            if ( state == true )
            {
                short pattern = 37449;
                glLineStipple( 1, pattern );
                glEnable( GL_LINE_STIPPLE );
            }

            if ( m_smooth )
            {
                glEnable(GL_BLEND);
                glEnable(GL_LINE_SMOOTH);
                glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
            }


        openvdb::Vec3d ptn, wmin, wmax;
        openvdb::Vec3s color;
        openvdb::CoordBBox bbox;
        //size_t cOffset = 0;

        openvdb::Vec3s sNodeColors[] = {
            openvdb::Vec3s(0.045, 0.045, 0.045),            // root
            openvdb::Vec3s(0.0432, 0.33, 0.0411023),        // first internal node level
            openvdb::Vec3s(0.871, 0.394, 0.01916),          // intermediate internal node levels
            openvdb::Vec3s(0.00608299, 0.279541, 0.625)     // leaf nodes
            };

        glBegin( GL_LINES );

        for (typename GridType::TreeType::NodeCIter iter = grid->tree().cbeginNode(); iter; ++iter)
        {
            iter.getBoundingBox(bbox);

            // Nodes are rendered as cell-centered
            const openvdb::Vec3d min(bbox.min().x()-0.5, bbox.min().y()-0.5, bbox.min().z()-0.5);
            const openvdb::Vec3d max(bbox.max().x()+0.5, bbox.max().y()+0.5, bbox.max().z()+0.5);

            //convert to world space
            wmin = grid->indexToWorld(min);
            wmax = grid->indexToWorld(max);

            //get colour based on level
            const int level = iter.getLevel();
            color = sNodeColors[(level == 0) ? 3 : (level == 1) ? 2 : 1];

            /*for (size_t n = 0; n < 8; ++n) {
                colors[cOffset++] = color[0];
                colors[cOffset++] = color[1];
                colors[cOffset++] = color[2];
            }*/

            glColor3f( color[0], color[1], color[2] );
           

            glVertex3f( wmin.x(), wmax.y(), wmax.z() );
            glVertex3f( wmax.x(), wmax.y(), wmax.z() );

            glVertex3f( wmin.x(), wmax.y(), wmin.z() );
            glVertex3f( wmax.x(), wmax.y(), wmin.z() );

            glVertex3f( wmin.x(), wmin.y(), wmax.z() );
            glVertex3f( wmax.x(), wmin.y(), wmax.z() );

            glVertex3f( wmin.x(), wmin.y(), wmin.z());
            glVertex3f( wmax.x(), wmin.y(), wmin.z() );

            glVertex3f( wmin.x(), wmin.y(), wmin.z() );
            glVertex3f( wmin.x(), wmin.y(), wmax.z() );

            glVertex3f( wmin.x(), wmax.y(), wmin.z() );
            glVertex3f( wmin.x(), wmax.y(), wmax.z() );

            glVertex3f( wmax.x(), wmin.y(), wmin.z() );
            glVertex3f( wmax.x(), wmin.y(), wmax.z() );

            glVertex3f( wmax.x(), wmax.y(), wmin.z() );
            glVertex3f( wmax.x(), wmax.y(), wmax.z() );

            glVertex3f( wmin.x(), wmin.y(), wmin.z() );
            glVertex3f( wmin.x(), wmax.y(), wmin.z() );

            glVertex3f( wmax.x(), wmin.y(), wmin.z() );
            glVertex3f( wmax.x(), wmax.y(), wmin.z() );

            glVertex3f( wmax.x(), wmin.y(), wmax.z() );
            glVertex3f( wmax.x(), wmax.y(), wmax.z() );

            glVertex3f( wmin.x(), wmin.y(), wmax.z() );
            glVertex3f( wmin.x(), wmax.y(), wmax.z() );

           
        } // end node iteration

            glEnd();
            glDisable( GL_LINE_STIPPLE );

            glPopAttrib();

        glEndList();

    }
    else
    {
        glCallList( m_treeID );
    }

}

template<typename GridType>
void vdbViewerNode::draw_meshDummy(typename GridType::Ptr grid)
{

}

template<typename GridType>
void vdbViewerNode::draw_mesh(typename GridType::Ptr grid)
{
    if ( m_meshID < 0 )
    {
        bool state = true;
        m_meshID = glGenLists(1);
        glNewList( m_meshID, GL_COMPILE_AND_EXECUTE );

            glPushAttrib( GL_ALL_ATTRIB_BITS);

            if ( state == true )
            {
                short pattern = 37449;
                glLineStipple( 1, pattern );
                glEnable( GL_LINE_STIPPLE );
            }

            if ( m_smooth )
            {
                glEnable(GL_BLEND);
                glEnable(GL_LINE_SMOOTH);
                glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
            }


        openvdb::FloatGrid::ConstPtr scalarGrid = openvdb::gridPtrCast<openvdb::FloatGrid>(grid);

        openvdb::tools::VolumeToMesh mesher(0.0);
        /*mesher(&scalarGrid);

        
        // Copy points and generate point normals.
        std::vector<GLfloat> points, normals(mesher.pointListSize() * 3);
        points.reserve(mesher.pointListSize() * 3);
        //normals.reserve(mesher.pointListSize() * 3);

        openvdb::tree::ValueAccessor<const typename GridType::TreeType> acc(grid->tree());
        typedef openvdb::math::Gradient<openvdb::math::GenericMap, openvdb::math::CD_2ND> Gradient;
        openvdb::math::GenericMap map(grid->transform());
        openvdb::Coord ijk;

        for (size_t n = 0, N = mesher.pointListSize(); n < N; ++n) {
            const openvdb::Vec3s& p = mesher.pointList()[n];
            points.push_back(p[0]);
            points.push_back(p[1]);
            points.push_back(p[2]);
        }

        // Copy primitives
        openvdb::tools::PolygonPoolList& polygonPoolList = mesher.polygonPoolList();
        size_t numQuads = 0;
        for (size_t n = 0, N = mesher.polygonPoolListSize(); n < N; ++n) {
            numQuads += polygonPoolList[n].numQuads();
        }

        std::vector<GLuint> indices;
        indices.reserve(numQuads * 4);
        openvdb::Vec3d normal, e1, e2;

        for (size_t n = 0, N = mesher.polygonPoolListSize(); n < N; ++n) {
            const openvdb::tools::PolygonPool& polygons = polygonPoolList[n];
            for (size_t i = 0, I = polygons.numQuads(); i < I; ++i) {
                const openvdb::Vec4I& quad = polygons.quad(i);
                indices.push_back(quad[0]);
                indices.push_back(quad[1]);
                indices.push_back(quad[2]);
                indices.push_back(quad[3]);

                e1 = mesher.pointList()[quad[1]];
                e1 -= mesher.pointList()[quad[0]];
                e2 = mesher.pointList()[quad[2]];
                e2 -= mesher.pointList()[quad[1]];
                normal = e1.cross(e2);

                const double length = normal.length();
                if (length > 1.0e-7) normal *= (1.0 / length);

                for (size_t v = 0; v < 4; ++v) {
                    normals[quad[v]*3]    = -normal[0];
                    normals[quad[v]*3+1]  = -normal[1];
                    normals[quad[v]*3+2]  = -normal[2];
                }
            }
        }


              glEnd();
            glDisable( GL_LINE_STIPPLE );

            glPopAttrib();

        glEndList()
 */          

    }
    else
    {
        glCallList( m_meshID );
    }

}

//
//  vdbViewerNode::draw_cropBox
//
void vdbViewerNode::draw_cropBox( bool state )
{
	if ( m_cropboxID < 0 )
	{

		m_cropboxID = glGenLists(1);
		glNewList( m_cropboxID, GL_COMPILE_AND_EXECUTE );

			glPushAttrib( GL_ALL_ATTRIB_BITS);

				if ( state == true )
				{
					short pattern = 37449;
					glLineStipple( 1, pattern );
					glEnable( GL_LINE_STIPPLE );
				}

				if ( m_smooth )
				{
					glEnable(GL_BLEND);
					glEnable(GL_LINE_SMOOTH);
					glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
				}

				glBegin( GL_LINES );

					glVertex3f( m_cropBoxLocal.min().x, m_cropBoxLocal.max().y, m_cropBoxLocal.max().z );
					glVertex3f( m_cropBoxLocal.max().x, m_cropBoxLocal.max().y, m_cropBoxLocal.max().z );

					glVertex3f( m_cropBoxLocal.min().x, m_cropBoxLocal.max().y, m_cropBoxLocal.min().z );
					glVertex3f( m_cropBoxLocal.max().x, m_cropBoxLocal.max().y, m_cropBoxLocal.min().z );

					glVertex3f( m_cropBoxLocal.min().x, m_cropBoxLocal.min().y, m_cropBoxLocal.max().z );
					glVertex3f( m_cropBoxLocal.max().x, m_cropBoxLocal.min().y, m_cropBoxLocal.max().z );

					glVertex3f( m_cropBoxLocal.min().x, m_cropBoxLocal.min().y, m_cropBoxLocal.min().z );
					glVertex3f( m_cropBoxLocal.max().x, m_cropBoxLocal.min().y, m_cropBoxLocal.min().z );

					glVertex3f( m_cropBoxLocal.min().x, m_cropBoxLocal.min().y, m_cropBoxLocal.min().z );
					glVertex3f( m_cropBoxLocal.min().x, m_cropBoxLocal.min().y, m_cropBoxLocal.max().z );

					glVertex3f( m_cropBoxLocal.min().x, m_cropBoxLocal.max().y, m_cropBoxLocal.min().z );
					glVertex3f( m_cropBoxLocal.min().x, m_cropBoxLocal.max().y, m_cropBoxLocal.max().z );

					glVertex3f( m_cropBoxLocal.max().x, m_cropBoxLocal.min().y, m_cropBoxLocal.min().z );
					glVertex3f( m_cropBoxLocal.max().x, m_cropBoxLocal.min().y, m_cropBoxLocal.max().z );

					glVertex3f( m_cropBoxLocal.max().x, m_cropBoxLocal.max().y, m_cropBoxLocal.min().z );
					glVertex3f( m_cropBoxLocal.max().x, m_cropBoxLocal.max().y, m_cropBoxLocal.max().z );

					glVertex3f( m_cropBoxLocal.min().x, m_cropBoxLocal.min().y, m_cropBoxLocal.min().z );
					glVertex3f( m_cropBoxLocal.min().x, m_cropBoxLocal.max().y, m_cropBoxLocal.min().z );

					glVertex3f( m_cropBoxLocal.max().x, m_cropBoxLocal.min().y, m_cropBoxLocal.min().z );
					glVertex3f( m_cropBoxLocal.max().x, m_cropBoxLocal.max().y, m_cropBoxLocal.min().z );

					glVertex3f( m_cropBoxLocal.max().x, m_cropBoxLocal.min().y, m_cropBoxLocal.max().z );
					glVertex3f( m_cropBoxLocal.max().x, m_cropBoxLocal.max().y, m_cropBoxLocal.max().z );

					glVertex3f( m_cropBoxLocal.min().x, m_cropBoxLocal.min().y, m_cropBoxLocal.max().z );
					glVertex3f( m_cropBoxLocal.min().x, m_cropBoxLocal.max().y, m_cropBoxLocal.max().z );

				glEnd();

				glDisable( GL_LINE_STIPPLE );
				glLineWidth(3.0);

				float w = m_cropBoxLocal.width()  * 0.2;
				float h = m_cropBoxLocal.height() * 0.2;
				float d = m_cropBoxLocal.depth()  * 0.2;

				// min coord sys
				glBegin( GL_LINES );
					glColor3f( 1.0, 0.0, 0.0 );
					glVertex3f( m_cropBoxLocal.min().x, m_cropBoxLocal.min().y, m_cropBoxLocal.min().z );
					glVertex3f( m_cropBoxLocal.min().x+w, m_cropBoxLocal.min().y, m_cropBoxLocal.min().z );

					glColor3f( 0.0, 1.0, 0.0 );
					glVertex3f( m_cropBoxLocal.min().x, m_cropBoxLocal.min().y, m_cropBoxLocal.min().z );
					glVertex3f( m_cropBoxLocal.min().x, m_cropBoxLocal.min().y+h, m_cropBoxLocal.min().z );

					glColor3f( 0.0, 0.0, 1.0 );
					glVertex3f( m_cropBoxLocal.min().x, m_cropBoxLocal.min().y, m_cropBoxLocal.min().z );
					glVertex3f( m_cropBoxLocal.min().x, m_cropBoxLocal.min().y, m_cropBoxLocal.min().z+d );
				glEnd();

				// max coord sys
				glBegin( GL_LINES );
					glColor3f( 1.0, 0.0, 0.0 );
					glVertex3f( m_cropBoxLocal.max().x, m_cropBoxLocal.max().y, m_cropBoxLocal.max().z );
					glVertex3f( m_cropBoxLocal.max().x-w, m_cropBoxLocal.max().y, m_cropBoxLocal.max().z );

					glColor3f( 0.0, 1.0, 0.0 );
					glVertex3f( m_cropBoxLocal.max().x, m_cropBoxLocal.max().y, m_cropBoxLocal.max().z );
					glVertex3f( m_cropBoxLocal.max().x, m_cropBoxLocal.max().y-h, m_cropBoxLocal.max().z );

					glColor3f( 0.0, 0.0, 1.0 );
					glVertex3f( m_cropBoxLocal.max().x, m_cropBoxLocal.max().y, m_cropBoxLocal.max().z );
					glVertex3f( m_cropBoxLocal.max().x, m_cropBoxLocal.max().y, m_cropBoxLocal.max().z-d );
				glEnd();

			glPopAttrib();
		glEndList();

	}
	else
	{
		glCallList( m_cropboxID );
	}

}


inline bool vdbViewerNode::filterPoint( MVector &val, float &radius, bool load )
{

	if ( m_filterMode == FILTER_PASSTHROUGH ||
		( load && m_filterOn == FILTER_ON_DRAW )
		) return true;

	switch( m_filterVariable )
	{
		case VALUE_FILTER:
		{
			switch( m_filterMode )
			{
				case FILTER_LESSTHAN:
				{
					return ( val.x < m_filterValc1.x &&
							 val.y < m_filterValc1.y &&
							 val.z < m_filterValc1.z );
				}
				case FILTER_MORETHAN:
				{
					return ( val.x > m_filterValc1.x &&
							 val.y > m_filterValc1.y &&
							 val.z > m_filterValc1.z );
				}
				case FILTER_EQUAL:
				{
					return ( val.x == m_filterValc1.x &&
							 val.y == m_filterValc1.y &&
							 val.z == m_filterValc1.z );
				}
				case FILTER_BETWEEN:
				{
					return ( ( val.x > m_filterValc1.x && val.x < m_filterValc2.x ) &&
							 ( val.y > m_filterValc1.y && val.y < m_filterValc2.y ) &&
							 ( val.z > m_filterValc1.z && val.z < m_filterValc2.z )  );
				}
				ERRV("Unknown value filter mode: ",m_filterMode);
			};
			break;
		}
		case RADIUS_FILTER:
		{
			switch( m_filterMode )
			{
				case FILTER_LESSTHAN:
				{
					if ( radius < m_filterValf1 ) return true;
					else return false;
				}
				case FILTER_MORETHAN:
				{
					if ( radius > m_filterValf1 ) return true;
					else return false;
				}
				case FILTER_EQUAL:
				{
					if ( radius == m_filterValf1 ) return true;
					else return false;
				}
				case FILTER_BETWEEN:
				{
					if ( radius > m_filterValf1 && radius < m_filterValf2 ) return true;
					else return false;
				}
				ERRV("Unknown filter mode: ",m_filterMode);
			};
			break;
		}
	};

	ERRV("Unknown filter variable: ",m_filterVariable);
	return false;
}


//
//  vdbViewerNode::draw
//
void vdbViewerNode::draw(   M3dView & view, const MDagPath & path,
							  M3dView::DisplayStyle style,
							  M3dView::DisplayStatus displaystatus )
{
	MStatus status;
	if ( !m_bDummyConnected || !m_bFilterChangeConnected ) connectDummyToShear();
	if ( !m_timeConnected ) connectToTime();


	MObject thisNode = thisMObject();
	MFnDependencyNode nodeFn(thisNode);

	// Start drawing
	view.beginGL();
		glPushMatrix();
		glPushAttrib( GL_ALL_ATTRIB_BITS);


		// display the bounding box when selected.
		if ( displaystatus != M3dView::kDormant )
		{

			if ( m_doCropBox )
			{
				draw_cropBox( false );
				draw_bbox( true );
			}
			else
			{
				draw_cropBox( true );
				draw_bbox( false );
			}

		}

        

        //show VDB Treee  
        if(m_showVoxelTree)
        {

            for (openvdb::GridPtrVec::const_iterator it = grids->begin(); it != grids->end(); ++it) {
                openvdb::GridBase::Ptr curGrid = *it;
                if (!curGrid) continue;
                if(curGrid->getName()==m_varnames[m_chan])
                {
                  processTypedGridTree(curGrid);
                }

            }
        }




        //TODO - write a helper that draws the VDB mesh for (float grid) signed distance fields 
        /*if(m_showmesh)
        {
            for (openvdb::GridPtrVec::const_iterator it = grids->begin(); it != grids->end(); ++it) {
                openvdb::GridBase::Ptr curGrid = *it;
                if (!curGrid) continue;

                if(curGrid->getName()==m_varnames[m_chan])
                {
                  processTypedGridMesh(curGrid);
                }
            }
        }*/


		bool vectorData = m_channelType == "vector" || m_channelType == "Vector";

		if ( m_pcID < 0 || (m_viewAs == 1) || (m_viewAs == 2) )
		{
			int slices = 8;
			if ( m_viewAs == 1 )
			{
				slices = m_circleSlices;
			}
			else
			{
				slices = m_diskSlices;
			}


			{
				m_diskID = glGenLists(1);
				glNewList( m_diskID, GL_COMPILE );
					if ( m_viewAs == 1 ) {
						if ( m_smooth ) {
							glEnable(GL_BLEND);
							glEnable(GL_LINE_SMOOTH);
							glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
						}
						gluQuadricDrawStyle(m_quadricObj,GLU_SILHOUETTE);
					} else if ( m_viewAs == 2 )
						gluQuadricDrawStyle(m_quadricObj,GLU_FILL);
					gluQuadricNormals(m_quadricObj,GLU_NONE);
					gluDisk( m_quadricObj, 0, 1.0, slices, 1 );
				glEndList();
				DEBUG(" + Create disk");
			}

			m_pcID = glGenLists(1);

			glNewList( m_pcID, GL_COMPILE_AND_EXECUTE );

				if ( m_smooth )
				{
					glEnable(GL_DEPTH_TEST);
					glDepthFunc(GL_LESS);
					glEnable(GL_BLEND);
					glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
					glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
				}
				glEnable( GL_POINT_SMOOTH );


				glPointSize( m_ptSize );

                int len = m_nPointsLoaded;
				MVector col( 1.0, 1.0, 1.0 );

				if ( m_viewAs == 0 )
				{
					glBegin( GL_POINTS );


						if ( vectorData )
						{
							// Vector and VECTOR channels
							if ( m_invert ) glColor3f( 1.0, 1.0, 1.0 );
							else glColor3f( 0.5, 0.5, 0.5 );
							for ( int i=0; i<len; i++ )
							{
								glVertex3f( m_pts[i].x, m_pts[i].y, m_pts[i].z );
							}
						}
						else
						{
							// FLOAT and COLOR channels
							for ( int i=0; i<len; i++ )
							{
								col = m_ptsColor[i];
								if ( m_filterOn == FILTER_ON_LOAD ||
									filterPoint( col, m_ptsRadius[i], false ) )
								{
									if ( m_abs )        col = MVector( fabs(col.x), fabs(col.y), fabs(col.z) );
									if ( m_invert )     col = MVector( 1.0,1.0,1.0 ) - col;
									if ( m_exp != 1.0 ) col = m_exp * col;
									glColor3f( col.x, col.y, col.z );
									glVertex3f( m_pts[i].x, m_pts[i].y, m_pts[i].z );
								}
							}
						}
					glEnd();
				}
				else if ( m_viewAs == 1 || m_viewAs == 2 )
				{
					if ( vectorData )
					{
						// Vector and VECTOR channels
						// do nothing
					}
					else
					{
						// FLOAT and COLOR channels
						float m[16];
						for ( int i=0; i<len; i++ )
						{
							col = m_ptsColor[i];
							if ( m_filterOn == FILTER_ON_LOAD ||
								filterPoint( col, m_ptsRadius[i], false ) )
							{
								if ( m_abs )        col = MVector( fabs(col.x), fabs(col.y), fabs(col.z) );
								if ( m_invert )     col = MVector( 1.0,1.0,1.0 ) - col;
								if ( m_exp != 1.0 ) col = m_exp * col;

                                MVector n( m_ptsVector[i].x, m_ptsVector[i].y, m_ptsVector[i].z );

                                //if we can get camera orient disks to them
                                MDagPath activeCamera;
                                status = view.getCamera ( activeCamera  );
                                if ( status == MS::kSuccess )
                                {
                                MMatrix worldSpace = activeCamera.inclusiveMatrix(&status);
                                    if ( status == MS::kSuccess )
                                    {
                                        n = MVector( worldSpace[2][0], worldSpace[2][1], worldSpace[2][2] );
                                    }
                                    else
                                    {
                                       //std::cerr<<"could not get camera matrix\n" ;
                                    }
                                
                                }
                                else
                                {
                                    //std::cerr<<"could not get active camera dag path\n";
                                }




								n.normalize();
								MVector ref(1,0,0);
								MVector up = ref ^ n;
								up.normalize();
								MVector right = up ^ n;
								right.normalize();

								m[0]  = right.x;
								m[1]  = right.y;
								m[2]  = right.z;
								m[3]  = 0.0;
								m[4]  = up.x;
								m[5]  = up.y;
								m[6]  = up.z;
								m[7]  = 0.0;
								m[8]  = n.x;
								m[9]  = n.y;
								m[10] = n.z;
								m[11] = 0.0;
								m[12] = m_pts[i].x;
								m[13] = m_pts[i].y;
								m[14] = m_pts[i].z;
								m[15] = 1.0;

								glPushMatrix();
									glMultMatrixf(m);
									glColor3f( col.x, col.y, col.z );
									glPushMatrix();
										glScalef(m_ptsRadius[i]*m_ptSize, m_ptsRadius[i]*m_ptSize, m_ptsRadius[i]*m_ptSize);
										glCallList(m_diskID);
									glPopMatrix();
								glPopMatrix();
							}
						}
					}
				}


                if ( m_viewAs == 4 )
                {
                    printShortListing(1);
                }

				if ( m_Vectors && !vectorData )
				{

					glBegin( GL_LINES );
						for ( int i=0; i<len; i++ )
						{
							if ( m_filterOn == FILTER_ON_LOAD ||
								filterPoint( m_ptsColor[i], m_ptsRadius[i], false ) )
							{
								if ( m_invertVectorsColor )
									glColor3f( -m_ptsVector[i].x, -m_ptsVector[i].y, -m_ptsVector[i].z );
								else
									glColor3f( m_ptsVector[i].x, m_ptsVector[i].y, m_ptsVector[i].z );

								glVertex3f( m_pts[i].x, m_pts[i].y, m_pts[i].z );
								glVertex3f( m_pts[i].x + ( m_ptsVector[i].x * m_VectorsSize ),
								m_pts[i].y + ( m_ptsVector[i].y * m_VectorsSize ),
								m_pts[i].z + ( m_ptsVector[i].z * m_VectorsSize ) );
							}
						}
					glEnd();
				}


				if ( vectorData )
				{
					glBegin( GL_LINES );
						for ( int i=0; i<len; i++ )
						{
							if ( m_invertVectorsColor )
								glColor3f( -m_ptsColor[i].x, -m_ptsColor[i].y, -m_ptsColor[i].z );
							else
								glColor3f( m_ptsColor[i].x, m_ptsColor[i].y, m_ptsColor[i].z );

							glVertex3f( m_pts[i].x, m_pts[i].y, m_pts[i].z );
							glVertex3f( m_pts[i].x + ( m_ptsColor[i].x * m_VectorsSize ),
										m_pts[i].y + ( m_ptsColor[i].y * m_VectorsSize ),
										m_pts[i].z + ( m_ptsColor[i].z * m_VectorsSize ) );
						}
					glEnd();
				}
			glEndList();

		}
		else
		{
			DEBUG(" + use display list");
			glCallList( m_pcID );
		}

		if ( m_showValue && !vectorData )
		{
            int len = m_nPointsLoaded;
			glColor3f(0.4,0.4,0.4);
			MString label;
			MPoint  lpos;
			for ( int i=0; i<len; i++ )
			{
				if ( m_filterOn == FILTER_ON_LOAD ||
					filterPoint( m_ptsColor[i], m_ptsRadius[i], false ) )
				{
					label = " ( ";
					label += m_ptsColor[i].x;
					label += ", ";
					label += m_ptsColor[i].y;
					label += ", ";
					label += m_ptsColor[i].z;
					label += " )";
					lpos.x  = m_pts[i].x;
					lpos.y  = m_pts[i].y;
					lpos.z  = m_pts[i].z;
					view.drawText( label, lpos, M3dView::kLeft );
				}
			}
		}

	  glPopAttrib();
	glPopMatrix();
  view.endGL();

}


//
//  vdbViewerNode::isBounded
//
bool vdbViewerNode::isBounded() const
{
	return true;
}


//
//  vdbViewerNode::boundingBox
//
MBoundingBox vdbViewerNode::boundingBox() const
{
	if ( m_doCropBox )
	{
		return m_cropBoxLocal;
	}
	return m_bbox;
}


//
//  vdbViewerNode::creator
//
void* vdbViewerNode::creator()
{
	return new vdbViewerNode();
}


//
//  vdbViewerNode::initialize
//
MStatus vdbViewerNode::initialize()
{
	MFnTypedAttribute	matAttr;
	MFnTypedAttribute	strAttr;
	MFnNumericAttribute	numAttr;
	MFnUnitAttribute	uAttr;
	MFnEnumAttribute	eAttr;
	MStatus				stat;

	aFilterChange = numAttr.create( "filterChange", "fic", MFnNumericData::kFloat, 0.0, &stat );
	CHECK_MSTATUS( numAttr.setStorable(false) );
	CHECK_MSTATUS( numAttr.setKeyable(false) );
	CHECK_MSTATUS( numAttr.setHidden(true) );
	CHECK_MSTATUS( addAttribute( aFilterChange ) );

	aDummy = numAttr.create( "dummy", "dum", MFnNumericData::kFloat, 0.0, &stat );
	CHECK_MSTATUS( numAttr.setStorable(false) );
	CHECK_MSTATUS( numAttr.setKeyable(false) );
	CHECK_MSTATUS( numAttr.setHidden(true) );
	CHECK_MSTATUS( addAttribute( aDummy ) );


	aVdbFile = strAttr.create( MString("vdbFile"), MString("vdb"), MFnData::kString, aVdbFile, &stat );
	CHECK_MSTATUS( strAttr.setStorable(true) );
	CHECK_MSTATUS( strAttr.setKeyable(false) );
	CHECK_MSTATUS( strAttr.setConnectable(true) );
	CHECK_MSTATUS( addAttribute( aVdbFile ) );

	aChannels = numAttr.create( "channel", "ch", MFnNumericData::kInt, 0.0, &stat );
	CHECK_MSTATUS( numAttr.setStorable(true) );
	CHECK_MSTATUS( numAttr.setKeyable(true) );
	CHECK_MSTATUS( numAttr.setConnectable(true) );
	CHECK_MSTATUS( addAttribute( aChannels ) );

	aViewAs = eAttr.create( "displayAs", "va", 0, &stat );
	CHECK_MSTATUS( eAttr.setStorable(true) );
	CHECK_MSTATUS( eAttr.setKeyable(true) );
	CHECK_MSTATUS( eAttr.setConnectable(true) );
	CHECK_MSTATUS( eAttr.addField("Points",			0) );
	CHECK_MSTATUS( eAttr.addField("Circle",	1) );
	CHECK_MSTATUS( eAttr.addField("Disk",	2) );
    CHECK_MSTATUS( eAttr.addField("None (just bound)",  3) );
    CHECK_MSTATUS( eAttr.addField("Print stats (and bound)",    4) );
	CHECK_MSTATUS( addAttribute( aViewAs ) );


	aPercent = numAttr.create( "percentLoaded", "pl", MFnNumericData::kFloat, 25.0, &stat );
	CHECK_MSTATUS( numAttr.setMin( 0.0f ) );
	CHECK_MSTATUS( numAttr.setMax(  100.0f ) );
	CHECK_MSTATUS( numAttr.setStorable(true) );
	CHECK_MSTATUS( numAttr.setKeyable(true) );
	CHECK_MSTATUS( numAttr.setConnectable(true) );
	CHECK_MSTATUS( addAttribute( aPercent ) );

	aExposure = numAttr.create( "exposure", "exp", MFnNumericData::kFloat, 0.0, &stat );
	CHECK_MSTATUS( numAttr.setMin( -100.0f ) );
	CHECK_MSTATUS( numAttr.setMax(  100.0f ) );
	CHECK_MSTATUS( numAttr.setSoftMin( -10.0f ) );
	CHECK_MSTATUS( numAttr.setSoftMax(  10.0f ) );
	CHECK_MSTATUS( numAttr.setStorable(true) );
	CHECK_MSTATUS( numAttr.setKeyable(true) );
	CHECK_MSTATUS( numAttr.setConnectable(true) );
	CHECK_MSTATUS( addAttribute( aExposure ) );

	aPointSize = numAttr.create( "pointSize", "ps", MFnNumericData::kFloat, 1.0, &stat );
	CHECK_MSTATUS( numAttr.setMin(    1.0f ) );
	CHECK_MSTATUS( numAttr.setMax(  100.0f ) );
	CHECK_MSTATUS( numAttr.setSoftMin(   0.001f ) );
	CHECK_MSTATUS( numAttr.setSoftMax(  50.0f ) );
	CHECK_MSTATUS( numAttr.setStorable(true) );
	CHECK_MSTATUS( numAttr.setKeyable(true) );
	CHECK_MSTATUS( numAttr.setConnectable(true) );
	CHECK_MSTATUS( addAttribute( aPointSize ) );

	aSmooth = numAttr.create( "smooth", "sm", MFnNumericData::kBoolean, false, &stat );
	CHECK_MSTATUS( numAttr.setStorable(true) );
	CHECK_MSTATUS( numAttr.setKeyable(true) );
	CHECK_MSTATUS( numAttr.setConnectable(true) );
	CHECK_MSTATUS( addAttribute( aSmooth ) );

	aInvert = numAttr.create( "invert", "i", MFnNumericData::kBoolean, false, &stat );
	CHECK_MSTATUS( numAttr.setStorable(true) );
	CHECK_MSTATUS( numAttr.setKeyable(true) );
	CHECK_MSTATUS( numAttr.setConnectable(true) );
	CHECK_MSTATUS( addAttribute( aInvert ) );

	aAbsoluteValue = numAttr.create( "absoluteValue", "abs", MFnNumericData::kBoolean, false, &stat );
	CHECK_MSTATUS( numAttr.setStorable(true) );
	CHECK_MSTATUS( numAttr.setKeyable(true) );
	CHECK_MSTATUS( numAttr.setConnectable(true) );
	CHECK_MSTATUS( addAttribute( aAbsoluteValue ) );

    ashowVoxelTree = numAttr.create( "showVoxelTree", "svt", MFnNumericData::kBoolean, false, &stat );
    CHECK_MSTATUS( numAttr.setStorable(true) );
    CHECK_MSTATUS( numAttr.setKeyable(true) );
    CHECK_MSTATUS( numAttr.setConnectable(true) );
    CHECK_MSTATUS( addAttribute( ashowVoxelTree ) );

	aShowVectors = numAttr.create( "showVectors", "sn", MFnNumericData::kBoolean, false, &stat );
	CHECK_MSTATUS( numAttr.setStorable(true) );
	CHECK_MSTATUS( numAttr.setKeyable(true) );
	CHECK_MSTATUS( numAttr.setConnectable(true) );
	CHECK_MSTATUS( addAttribute( aShowVectors ) );

	aVectorsSize = numAttr.create( "VectorsSize", "ns", MFnNumericData::kFloat, 1.0, &stat );
	CHECK_MSTATUS( numAttr.setMin(  0.001f ) );
	CHECK_MSTATUS( numAttr.setMax(  10.0f ) );
	CHECK_MSTATUS( numAttr.setStorable(true) );
	CHECK_MSTATUS( numAttr.setKeyable(true) );
	CHECK_MSTATUS( numAttr.setConnectable(true) );
	CHECK_MSTATUS( addAttribute( aVectorsSize ) );

	aInvertVectorsColor = numAttr.create( "negateVectorsColor", "nnc", MFnNumericData::kBoolean, false, &stat );
	CHECK_MSTATUS( numAttr.setStorable(true) );
	CHECK_MSTATUS( numAttr.setKeyable(true) );
	CHECK_MSTATUS( numAttr.setConnectable(true) );
	CHECK_MSTATUS( addAttribute( aInvertVectorsColor ) );

	aUseCropBox = numAttr.create( "useCropBox", "cb", MFnNumericData::kBoolean, false, &stat );
	CHECK_MSTATUS( numAttr.setStorable(true) );
	CHECK_MSTATUS( numAttr.setKeyable(true) );
	CHECK_MSTATUS( numAttr.setConnectable(true) );
	CHECK_MSTATUS( addAttribute( aUseCropBox ) );

	aCropBoxMin = numAttr.createPoint( "cropBoxMin", "cmn", &stat );
	CHECK_MSTATUS( numAttr.setDefault( 0.5, 0.0, 0.0 ) );
	CHECK_MSTATUS( numAttr.setMin( 0.0, 0.0, 0.0 ) );
	CHECK_MSTATUS( numAttr.setMax( 1.0, 1.0, 1.0 ) );
	CHECK_MSTATUS( numAttr.setStorable(    true    ) );
	CHECK_MSTATUS( numAttr.setKeyable(     true    ) );
	CHECK_MSTATUS( numAttr.setConnectable( true    ) );
	CHECK_MSTATUS( addAttribute( aCropBoxMin ) );

	aCropBoxMax = numAttr.createPoint( "cropBoxMax", "cmx", &stat );
	CHECK_MSTATUS( numAttr.setDefault( 0.51, 1.0, 1.0 ) );
	CHECK_MSTATUS( numAttr.setMin( 0.0, 0.0, 0.0 ) );
	CHECK_MSTATUS( numAttr.setMax( 1.0, 1.0, 1.0 ) );
	CHECK_MSTATUS( numAttr.setStorable(    true    ) );
	CHECK_MSTATUS( numAttr.setKeyable(     true    ) );
	CHECK_MSTATUS( numAttr.setConnectable( true    ) );
	CHECK_MSTATUS( addAttribute( aCropBoxMax ) );

	aShowProgress = numAttr.create( "showProgressWindow", "spw", MFnNumericData::kBoolean, false, &stat );
	CHECK_MSTATUS( numAttr.setStorable(true) );
	CHECK_MSTATUS( numAttr.setKeyable(true) );
	CHECK_MSTATUS( numAttr.setConnectable(true) );
	CHECK_MSTATUS( addAttribute( aShowProgress ) );

	aTime = uAttr.create( "time", "tm", MFnUnitAttribute::kTime, 0.0 );
	uAttr.setReadable(true);
	uAttr.setWritable(true);
	uAttr.setStorable(true);
	uAttr.setKeyable(true);
	uAttr.setConnectable(true);
	stat = addAttribute( aTime );

	aShowValue = numAttr.create( "showValue", "sv", MFnNumericData::kBoolean, false, &stat );
	CHECK_MSTATUS( numAttr.setStorable(true) );
	CHECK_MSTATUS( numAttr.setKeyable(true) );
	CHECK_MSTATUS( numAttr.setConnectable(true) );
	CHECK_MSTATUS( addAttribute( aShowValue ) );

	aFilterVariable = eAttr.create( "filterVariable", "fv", 0, &stat );
	CHECK_MSTATUS( eAttr.setStorable(true) );
	CHECK_MSTATUS( eAttr.setKeyable(true) );
	CHECK_MSTATUS( eAttr.setConnectable(true) );
	CHECK_MSTATUS( eAttr.addField("Channel",	0) );
	CHECK_MSTATUS( eAttr.addField("Radius",		1) );
	CHECK_MSTATUS( addAttribute( aFilterVariable ) );

	aFilterMode = eAttr.create( "filterMode", "fm", 0, &stat );
	CHECK_MSTATUS( eAttr.setStorable(true) );
	CHECK_MSTATUS( eAttr.setKeyable(true) );
	CHECK_MSTATUS( eAttr.setConnectable(true) );
	CHECK_MSTATUS( eAttr.addField("Off",			0) );
	CHECK_MSTATUS( eAttr.addField("v < v1",			1) );
	CHECK_MSTATUS( eAttr.addField("v > v1",			2) );
	CHECK_MSTATUS( eAttr.addField("v == v1",		3) );
	CHECK_MSTATUS( eAttr.addField("v1 < v < v2",	4) );
	CHECK_MSTATUS( addAttribute( aFilterMode ) );

	aFilterValf1 = numAttr.create( "filterValf1", "fvf1", MFnNumericData::kFloat, 0.0, &stat );
	CHECK_MSTATUS( numAttr.setStorable(true) );
	CHECK_MSTATUS( numAttr.setKeyable(true) );
	CHECK_MSTATUS( numAttr.setConnectable(true) );
	CHECK_MSTATUS( addAttribute( aFilterValf1 ) );

	aFilterValf2 = numAttr.create( "filterValf2", "fvf2", MFnNumericData::kFloat, 1.0, &stat );
	CHECK_MSTATUS( numAttr.setStorable(true) );
	CHECK_MSTATUS( numAttr.setKeyable(true) );
	CHECK_MSTATUS( numAttr.setConnectable(true) );
	CHECK_MSTATUS( addAttribute( aFilterValf2 ) );

	aFilterValc1 = numAttr.createColor( "filterValc1", "fvc1", &stat );
	CHECK_MSTATUS( numAttr.setDefault(0.0, 0.0, 0.0) );
	CHECK_MSTATUS( numAttr.setStorable(true) );
	CHECK_MSTATUS( numAttr.setKeyable(true) );
	CHECK_MSTATUS( numAttr.setConnectable(true) );
	CHECK_MSTATUS( addAttribute( aFilterValc1 ) );

	aFilterValc2 = numAttr.createColor( "filterValc2", "fvc2", &stat );
	CHECK_MSTATUS( numAttr.setDefault(1.0, 1.0, 1.0) );
	CHECK_MSTATUS( numAttr.setStorable(true) );
	CHECK_MSTATUS( numAttr.setKeyable(true) );
	CHECK_MSTATUS( numAttr.setConnectable(true) );
	CHECK_MSTATUS( addAttribute( aFilterValc2 ) );

	aCircleSlices = numAttr.create( "circleSlices", "csl", MFnNumericData::kInt, 8, &stat );
	CHECK_MSTATUS( numAttr.setMin(3) );
	CHECK_MSTATUS( numAttr.setMax(32) );
	CHECK_MSTATUS( numAttr.setStorable(true) );
	CHECK_MSTATUS( numAttr.setKeyable(true) );
	CHECK_MSTATUS( numAttr.setConnectable(true) );
	CHECK_MSTATUS( addAttribute( aCircleSlices ) );

	aDiskSlices = numAttr.create( "diskSlices", "dsl", MFnNumericData::kInt, 8, &stat );
	CHECK_MSTATUS( numAttr.setMin(3) );
	CHECK_MSTATUS( numAttr.setMax(32) );
	CHECK_MSTATUS( numAttr.setStorable(true) );
	CHECK_MSTATUS( numAttr.setKeyable(true) );
	CHECK_MSTATUS( numAttr.setConnectable(true) );
	CHECK_MSTATUS( addAttribute( aDiskSlices ) );

	aFilterOn = eAttr.create( "filterOn", "fo", 1, &stat );
	CHECK_MSTATUS( eAttr.setStorable(true) );
	CHECK_MSTATUS( eAttr.setKeyable(true) );
	CHECK_MSTATUS( eAttr.setConnectable(true) );
	CHECK_MSTATUS( eAttr.addField("Load",			0) );
	CHECK_MSTATUS( eAttr.addField("Draw",			1) );
	CHECK_MSTATUS( addAttribute( aFilterOn ) );

	// output attr
	aChannelNames = strAttr.create(
			MString("channelNames"),
			MString("chn"),
			MFnData::kStringArray,
			aChannelNames,
			&stat
			);
	CHECK_MSTATUS( strAttr.setHidden(      true  ) );
	CHECK_MSTATUS( strAttr.setStorable(    false ) );
	CHECK_MSTATUS( strAttr.setKeyable(     false ) );
	CHECK_MSTATUS( strAttr.setConnectable( true  ) );
	CHECK_MSTATUS( addAttribute( aChannelNames ) );

	aNumPoints = numAttr.create( "numPoints", "np", MFnNumericData::kInt, 0.0, &stat );
	CHECK_MSTATUS( numAttr.setHidden(      true  ) );
	CHECK_MSTATUS( numAttr.setStorable(    false ) );
	CHECK_MSTATUS( numAttr.setKeyable(     false ) );
	CHECK_MSTATUS( numAttr.setConnectable( true  ) );
	CHECK_MSTATUS( addAttribute( aNumPoints ) );

	aNumPointsLoaded = numAttr.create( "numPointsLoaded", "npl", MFnNumericData::kInt, 0.0, &stat );
	CHECK_MSTATUS( numAttr.setHidden(      true  ) );
	CHECK_MSTATUS( numAttr.setStorable(    false ) );
	CHECK_MSTATUS( numAttr.setKeyable(     false ) );
	CHECK_MSTATUS( numAttr.setConnectable( true  ) );
	CHECK_MSTATUS( addAttribute( aNumPointsLoaded ) );


	// relationships
	//
	CHECK_MSTATUS( attributeAffects( aFilterVariable, aFilterChange ) );
	CHECK_MSTATUS( attributeAffects( aFilterMode,     aFilterChange ) );
	CHECK_MSTATUS( attributeAffects( aFilterValf1,    aFilterChange ) );
	CHECK_MSTATUS( attributeAffects( aFilterValf2,    aFilterChange ) );
	CHECK_MSTATUS( attributeAffects( aFilterValc1,    aFilterChange ) );
	CHECK_MSTATUS( attributeAffects( aFilterValc2,    aFilterChange ) );
	CHECK_MSTATUS( attributeAffects( aFilterOn,       aFilterChange ) );

	CHECK_MSTATUS( attributeAffects( aVdbFile,            aDummy ) );
	CHECK_MSTATUS( attributeAffects( aChannels,           aDummy ) );
	CHECK_MSTATUS( attributeAffects( aViewAs,             aDummy ) );
	CHECK_MSTATUS( attributeAffects( aExposure,           aDummy ) );
	CHECK_MSTATUS( attributeAffects( aPointSize,          aDummy ) );
	CHECK_MSTATUS( attributeAffects( aSmooth,             aDummy ) );
	CHECK_MSTATUS( attributeAffects( aInvert,             aDummy ) );
	CHECK_MSTATUS( attributeAffects( aAbsoluteValue,      aDummy ) );
	CHECK_MSTATUS( attributeAffects( aShowVectors,        aDummy ) );
    CHECK_MSTATUS( attributeAffects( ashowVoxelTree,      aDummy ) );
	CHECK_MSTATUS( attributeAffects( aVectorsSize,        aDummy ) );
	CHECK_MSTATUS( attributeAffects( aInvertVectorsColor, aDummy ) );
	CHECK_MSTATUS( attributeAffects( aPercent,            aDummy ) );
	CHECK_MSTATUS( attributeAffects( aUseCropBox,         aDummy ) );
	CHECK_MSTATUS( attributeAffects( aCropBoxMin,         aDummy ) );
	CHECK_MSTATUS( attributeAffects( aCropBoxMax,         aDummy ) );
	CHECK_MSTATUS( attributeAffects( aShowProgress,       aDummy ) );
	CHECK_MSTATUS( attributeAffects( aTime,               aDummy ) );
	CHECK_MSTATUS( attributeAffects( aShowValue,          aDummy ) );
	CHECK_MSTATUS( attributeAffects( aCircleSlices,       aDummy ) );
	CHECK_MSTATUS( attributeAffects( aDiskSlices,         aDummy ) );

	return MS::kSuccess;
}


//
//  vdbViewerNode::connectDummyToShear
//
MStatus vdbViewerNode::connectDummyToShear()
{
	MStatus status;

	MFnDagNode dagFn(thisMObject(), &status);
	CHECKERR(status,"dagFn");
	MObject parent = dagFn.parent(0, &status);
	if(parent.isNull())
		return MS::kFailure;

	dagFn.setObject(parent);
	MPlug shearXYPlug = dagFn.findPlug("shearXY",&status);
	CHECKERR(status,"dagFn.findPlug(shearXY)");
	MPlug dummyPlug = MPlug(thisMObject(), aDummy);
	if(!dummyPlug.isConnected())
	{
		MDGModifier mod;
		mod.connect( dummyPlug, shearXYPlug );
		mod.doIt();
		m_bDummyConnected = true;
	}
	MPlug shearXZPlug = dagFn.findPlug("shearXZ",&status);
	CHECKERR(status,"dagFn.findPlug(shearXZ)");
	MPlug filterChangePlug = MPlug(thisMObject(), aFilterChange);
	if(!filterChangePlug.isConnected())
	{
		MDGModifier mod;
		mod.connect( filterChangePlug, shearXZPlug );
		mod.doIt();
		m_bFilterChangeConnected = true;
	}
	return status;

}


//
//  vdbViewerNode::connectToTime
//
MStatus vdbViewerNode::connectToTime()
{
	MStatus status;

	// our time plug
	MPlug timePlug( thisMObject(), aTime );

	MSelectionList sList;
	MGlobal::getSelectionListByName("time1", sList);
	unsigned int nMatches = sList.length();
	if ( nMatches > 0 )
	{
		MObject timeDepObj;
		sList.getDependNode(0, timeDepObj);
		MFnDependencyNode timeDep( timeDepObj );
		MPlug outTimePlug = timeDep.findPlug("outTime",&status);
		CHECKERR(status,"timeDep.findPlug(outTime)");
		if ( !timePlug.isConnected() )
		{
			MDGModifier mod;
			mod.connect( outTimePlug, timePlug );
			mod.doIt();
			m_timeConnected = true;
		}
		else
		{
			m_timeConnected = true;
		}
	}

	return status;

}





//
//  vdbViewerNode::loadVDB
//
MStatus vdbViewerNode::loadVDB()
{
	MStatus stat;
	float bbox[6];
    

    m_vdbFileName = buildFileName().asChar();


	openvdb::initialize();

    /// @todo Remove the following at some point:
    openvdb::Grid<openvdb::tree::Tree4<bool, 4, 3, 3>::Type>::registerGrid();
    openvdb::Grid<openvdb::tree::Tree4<float, 4, 3, 3>::Type>::registerGrid();
    openvdb::Grid<openvdb::tree::Tree4<double, 4, 3, 3>::Type>::registerGrid();
    openvdb::Grid<openvdb::tree::Tree4<int32_t, 4, 3, 3>::Type>::registerGrid();
    openvdb::Grid<openvdb::tree::Tree4<int64_t, 4, 3, 3>::Type>::registerGrid();
    openvdb::Grid<openvdb::tree::Tree4<openvdb::Vec2i, 4, 3, 3>::Type>::registerGrid();
    openvdb::Grid<openvdb::tree::Tree4<openvdb::Vec2s, 4, 3, 3>::Type>::registerGrid();
    openvdb::Grid<openvdb::tree::Tree4<openvdb::Vec2d, 4, 3, 3>::Type>::registerGrid();
    openvdb::Grid<openvdb::tree::Tree4<openvdb::Vec3i, 4, 3, 3>::Type>::registerGrid();
    openvdb::Grid<openvdb::tree::Tree4<openvdb::Vec3f, 4, 3, 3>::Type>::registerGrid();
    openvdb::Grid<openvdb::tree::Tree4<openvdb::Vec3d, 4, 3, 3>::Type>::registerGrid();


	

    vdbOpenInitilise(m_vdbFileName,stat);
    if ( stat != MS::kSuccess )
        {
            MString err = "unable to open vdbfile: ";
            err += m_vdbFileName.c_str();
            ERR(err);
            return MS::kFailure;
        }
	


	#ifdef _DEBUG
	cout <<endl<<endl<<"loadVDB : channel "<<m_chan<<" of "<<m_vdbFileName.c_str()<<endl;
	#endif



	int nChannels = 0;

	// -char **vartypes = NULL;
    std::vector<std::string> vartypes;

    // get thenames and types of grids in the vdb file
    m_varnames.clear();
    vdbGetNumberNamesTypes( nChannels,m_varnames,vartypes);

	// avoid reading a channel that does not exist.
	//
	if ( m_chan >= nChannels ) m_chan = nChannels-1;

	// set the bounding box
	int worldspace =1;
	vdbGetBound( worldspace, bbox,stat);
    if ( stat != MS::kSuccess )
    {
            MString err = "no grids in vdbfile: ";
            err += m_vdbFileName.c_str();
            ERR(err);
            return MS::kFailure;
    }

	m_bbox = MBoundingBox( MPoint( bbox[0], bbox[1], bbox[2] ), MPoint( bbox[3], bbox[4], bbox[5] ) );

	compute_cropBox_world();

	
	// set the channel names
	//
	#ifdef _DEBUG
	cout <<"VDB Volume has: "<<nChannels<<" channels"<<endl;
	#endif


	m_channelNames.clear();


    for ( int i=0; i<nChannels; i++ )
    {
        MString name = m_varnames[i].c_str();
        name += " ";
        name += vartypes[i].c_str();
        m_channelNames.append( name );
        DEBUG("  + "+name);

    }

    
	// find how many points are available
	//
	
    nPoints = 0;
    vdbGetNumberPoints( nPoints);

	m_nPoints = nPoints;
	m_nPointsLoaded = 0;

	

	#ifdef _DEBUG
	cout <<"VDB Point cloud contains "<<nPoints<<" points"<<endl;
	#endif



	// which channel should we load ?
	//
	MString vType( vartypes[m_chan].c_str() );
	// - int isVector = ( vType == "float" )? 0:1;

	// set the current channel type
	//
	m_channelType = vType;

	// compute channel offset
	//
	int chanOffset = 0;
	for ( int ch=0; ch<m_chan; ch++ )
	{
		if ( MString(vartypes[ch].c_str()) == "float" )
			chanOffset += 1;
		else
			chanOffset += 3;
	}
	#ifdef _DEBUG
	cout <<"channel offset for "<<m_varnames[m_chan]<<" is "<<chanOffset<<" ( isVector == "<<isVector<<" : "<<vartypes[m_chan]<<" )"<<endl;
	#endif

	// clear the arrays
	//
	m_pts.clear();
	m_ptsColor.clear();
	m_ptsVector.clear();
	m_ptsRadius.clear();

	#ifdef _DEBUG
	cout <<"m_cropBox = "<<m_cropBox.min()<<" "<<m_cropBox.max()<<endl;
	#endif

	// init progress window
    showProgress = false;
    progressStep = 100.0;

    if(m_ShowProgress)
    {
	if ( MGlobal::mayaState() ==  MGlobal::kInteractive )
		showProgress = MProgressWindow::reserve();
    }


	if ( showProgress )
	{
		MProgressWindow:: setInterruptable(false);
		MString pmsg = "Loading ";
		pmsg += m_nPoints;
		pmsg += " points...";
		MProgressWindow:: setProgressStatus(pmsg);
		MProgressWindow:: setTitle("vdbViewer");
		MProgressWindow::setProgressMin(0);
		MProgressWindow::setProgressMax(m_nPoints);
		MProgressWindow::startProgress();
		progressStep = floor(m_nPoints / 100);
	}


	

		// how many points are we going to load ?
		//
		float factor = m_percent * 0.01f;
        numPts = (int) fmax( 1.0f, floorf((float) nPoints * factor) );


		#ifdef _DEBUG
		cout <<"Loading "<<numPts<<" of "<<nPoints<<" points ( "<<m_percent<<"% )"<<endl;
		#endif

		// allocate memory
		//

		#ifdef _DEBUG
		cout <<"Allocating "<<numPts<<" points..."<<endl;
		#endif

		stat = m_pts.setLength(numPts);
		if ( stat != MS::kSuccess )
		{
			MString err = "Memory allocation for ";
			err += numPts;
			err += " points failed";
			ERR(err);
			return MS::kFailure;
		}
		stat = m_ptsColor.setLength(numPts);
		if ( stat != MS::kSuccess )
		{
			MString err = "Memory allocation for ";
			err += numPts;
			err += " colors failed";
			ERR(err);
			return MS::kFailure;
		}
		stat = m_ptsVector.setLength(numPts);
		if ( stat != MS::kSuccess )
		{
			MString err = "Memory allocation for ";
			err += numPts;
			err += " Vectors failed";
			ERR(err);
			return MS::kFailure;
		}
		stat = m_ptsRadius.setLength(numPts);
		if ( stat != MS::kSuccess )
		{
			MString err = "Memory allocation for ";
			err += numPts;
			err += " radius failed";
			ERR(err);
			return MS::kFailure;
		}

        s = floor( float(nPoints) / float(numPts) );

		MVector pColor;
		//bool test = m_filterOn == FILTER_ON_DRAW;

		#ifdef _DEBUG
		cout <<"reload all:"<<endl;
		cout <<"   condition is "<<test<<" ( "<<m_filterOn<<" )"<<endl;
		cout <<"   filterVar is "<<m_filterVariable<<endl;
		cout <<"   filterMode is "<<m_filterMode<<endl;
		cout <<"   v1 = "<<m_filterValf1<<endl;
		cout <<"   v2 = "<<m_filterValf2<<endl;
		#endif


        // lopp grids pic the requested one to draw
        for (openvdb::GridPtrVec::const_iterator it = grids->begin(); it != grids->end(); ++it) {
            openvdb::GridBase::Ptr curGrid = *it;
            if (!curGrid) continue;

            if(curGrid->getName()==m_varnames[m_chan])
            {
            processTypedGrid(curGrid);
            }
        }
		

		#ifdef _DEBUG
		cout <<"   Loaded "<<m_nPointsLoaded<<" points !"<<endl;
		#endif
	

	if ( showProgress )
	{
		MProgressWindow:: endProgress();
	}



	return MS::kSuccess;
}








//
//  initializePlugin
//
MStatus initializePlugin( MObject obj )
{


	MFnPlugin plugin( obj, "Milo Green", "1.0", "any");

	CHECK_MSTATUS(
			plugin.registerNode(
					"vdbViewerNode",
					vdbViewerNode::id,
					&vdbViewerNode::creator,
					&vdbViewerNode::initialize,
					MPxNode::kLocatorNode
					)
			);

	MGlobal::executeCommandOnIdle( "rehash;eval \"source AEvdbViewerNodeTemplate.mel\";", false );

	return MS::kSuccess;
}


//
//  uninitializePlugin
//
MStatus uninitializePlugin( MObject obj)
{
	MFnPlugin plugin( obj );

	CHECK_MSTATUS( plugin.deregisterNode( vdbViewerNode::id ) );

	return MS::kSuccess;
}

/// END ////