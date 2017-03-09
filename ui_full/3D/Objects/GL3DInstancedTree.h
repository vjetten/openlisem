#ifndef GL3DINSTANCEDTREE_H
#define GL3DINSTANCEDTREE_H



#include <3D/GL3DWidget.h>
#include <3D/Graphics/GL3DModels.h>
#include <3D/World/GL3DSurface.h>
#include <3D/Objects/GL3DGLInstancedModel.h>

class GL3DInstancedTree : public GL3DGLInstancedModel
{

public:
    GL3DInstancedTree() : GL3DGLInstancedModel()
    {
        //this->SetModelHighp("tree_lowp/tree_low.obj",100);
        this->SetSmoothing(false,0,0);
        this->SetIncrementDistance(20.0);
        this->SetRandomParameters(2.5,0.25);
        this->SetModelLowp("tree_lowp/tree_low.obj",2000);
        this->SetMaxInstances(10000,100);
    }

};


#endif // GL3DINSTANCEDTREE_H
