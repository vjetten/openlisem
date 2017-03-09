#ifndef GL3DGRASS_H
#define GL3DGRASS_H

#include <3D/GL3DWidget.h>
#include <3D/Graphics/GL3DModels.h>
#include <3D/World/GL3DSurface.h>
#include <3D/Objects/GL3DGLInstancedModel.h>

class GL3DGrass : public GL3DGLInstancedModel
{

public:
    GL3DGrass() : GL3DGLInstancedModel()
    {

        this->SetModelHighp("grass_highp/grass_low_poly.obj",200,0.2);
        this->SetSmoothing(false,0,0);
        this->SetIncrementDistance(2.0);
        this->SetRandomParameters(2.5,0.25);
        //this->SetModelLowp("grass_highp/grass_low_poly.obj",100,0.2);
        this->SetMaxInstances(2000,100);
    }

};


#endif // GL3DGRASS_H
