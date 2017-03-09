#include "3D/World/GL3DWorld.h"
#include "3D/World/GL3DCamera.h"
#include "3D/Graphics/GL3DShaders.h"
#include "3D/Graphics/GL3DTextures.h"
#include "3D/Graphics/GL3DGeometry.h"
#include "ui_full/3D/Objects/GL3DObject.h"
#include "3D/Graphics/GL3DMapMath.h"
#include "3D/Graphics/GL3DModels.h"
#include "3D/World/GL3DSurface.h"

#include "3D/Graphics/GL3DDrawFunctions.h"

 void GL3DDrawFunctions::BindGeometryInstanced(GL3DWidget * gl, QOpenGLVertexArrayObject &object,GL3DShader * s, GL3DGeometry * g, QOpenGLBuffer &MatrixBuffer)
{


    if(!g->is_2d_data)
    {
        if(g->uses_index)
        {
            object.bind();
            g->m_vertex.bind();
            g->m_index.bind();

            s->m_program->bind();

            int vertexlocation = s->m_program->attributeLocation("posAttr");
            int colorlocation = s->m_program->attributeLocation("colAttr");
            int uvlocation = s->m_program->attributeLocation("uvcAttr");
            int normallocation = s->m_program->attributeLocation("norAttr");
            int tangentlocation = s->m_program->attributeLocation("tanAttr");
            int bitangentlocation = s->m_program->attributeLocation("btaAttr");

            s->m_program->enableAttributeArray(vertexlocation);
            s->m_program->enableAttributeArray(colorlocation);

            s->m_program->enableAttributeArray(uvlocation);
            s->m_program->enableAttributeArray(normallocation);

            s->m_program->enableAttributeArray(tangentlocation);
            s->m_program->enableAttributeArray(bitangentlocation);

            s->m_program->setAttributeBuffer(vertexlocation, GL_FLOAT, Vertex::positionOffset(), Vertex::PositionTupleSize, Vertex::stride());
            s->m_program->setAttributeBuffer(colorlocation, GL_FLOAT, Vertex::colorOffset(), Vertex::ColorTupleSize, Vertex::stride());
            s->m_program->setAttributeBuffer(normallocation, GL_FLOAT, Vertex::normalOffset(), Vertex::NormalTupleSize, Vertex::stride());
            s->m_program->setAttributeBuffer(uvlocation, GL_FLOAT, Vertex::uvOffset(), Vertex::UVTupleSize, Vertex::stride());
            s->m_program->setAttributeBuffer(tangentlocation, GL_FLOAT, Vertex::tangentOffset(), Vertex::TangentTupleSize, Vertex::stride());
            s->m_program->setAttributeBuffer(bitangentlocation, GL_FLOAT, Vertex::bitangentOffset(), Vertex::BiTangentTupleSize, Vertex::stride());

            MatrixBuffer.bind();

            int matrixlocation = s->m_program->attributeLocation("matAttr");
            s->m_program->enableAttributeArray(matrixlocation+ 0);
            s->m_program->setAttributeBuffer(matrixlocation+ 0, GL_FLOAT, 0 * sizeof(QVector3D), 3, 1 * sizeof(QVector3D));
            gl->gl->glVertexAttribDivisor(matrixlocation + 0, 1);

            int e = gl->gl->glGetError();
            if(e != GL_NO_ERROR)
            {
                qDebug() << "error in binding vao 0 " << e;
            }

            object.release();
            MatrixBuffer.release();
            g->m_vertex.release();
            g->m_index.release();
            s->m_program->release();
        }
    }
}


 void GL3DDrawFunctions::BindGeometry(GL3DWidget * gl, QOpenGLVertexArrayObject &object,GL3DShader * s, GL3DGeometry * g)
{
    if(!g->is_2d_data)
    {
        if(g->uses_index)
        {
            object.bind();
            g->m_vertex.bind();
            g->m_index.bind();

            s->m_program->bind();

            int vertexlocation = s->m_program->attributeLocation("posAttr");
            int colorlocation = s->m_program->attributeLocation("colAttr");
            int uvlocation = s->m_program->attributeLocation("uvcAttr");
            int normallocation = s->m_program->attributeLocation("norAttr");
            int tangentlocation = s->m_program->attributeLocation("tanAttr");
            int bitangentlocation = s->m_program->attributeLocation("btaAttr");

            s->m_program->enableAttributeArray(vertexlocation);
            s->m_program->enableAttributeArray(colorlocation);

            s->m_program->enableAttributeArray(uvlocation);
            s->m_program->enableAttributeArray(normallocation);

            s->m_program->enableAttributeArray(tangentlocation);
            s->m_program->enableAttributeArray(bitangentlocation);

            s->m_program->setAttributeBuffer(vertexlocation, GL_FLOAT, Vertex::positionOffset(), Vertex::PositionTupleSize, Vertex::stride());
            s->m_program->setAttributeBuffer(colorlocation, GL_FLOAT, Vertex::colorOffset(), Vertex::ColorTupleSize, Vertex::stride());
            s->m_program->setAttributeBuffer(normallocation, GL_FLOAT, Vertex::normalOffset(), Vertex::NormalTupleSize, Vertex::stride());
            s->m_program->setAttributeBuffer(uvlocation, GL_FLOAT, Vertex::uvOffset(), Vertex::UVTupleSize, Vertex::stride());
            s->m_program->setAttributeBuffer(tangentlocation, GL_FLOAT, Vertex::tangentOffset(), Vertex::TangentTupleSize, Vertex::stride());
            s->m_program->setAttributeBuffer(bitangentlocation, GL_FLOAT, Vertex::bitangentOffset(), Vertex::BiTangentTupleSize, Vertex::stride());

            object.release();
            g->m_vertex.release();
            g->m_index.release();
            s->m_program->release();
        }else
        {

            object.bind();
            g->m_vertex.bind();

            s->m_program->bind();

            int vertexlocation = s->m_program->attributeLocation("posAttr");
            int colorlocation = s->m_program->attributeLocation("colAttr");
            int texturelocation = 0;

            s->m_program->enableAttributeArray(vertexlocation);
            s->m_program->enableAttributeArray(colorlocation);

            s->m_program->setAttributeBuffer(vertexlocation, GL_FLOAT, Vertex::positionOffset(), Vertex::PositionTupleSize, Vertex::stride());
            s->m_program->setAttributeBuffer(colorlocation, GL_FLOAT, Vertex::colorOffset(), Vertex::ColorTupleSize, Vertex::stride());

            object.release();
            g->m_vertex.release();
            s->m_program->release();

        }
    }else
    {
        object.bind();
        g->m_vertex.bind();

        s->m_program->bind();

        int vertexlocation = s->m_program->attributeLocation("posAttr");
        int texturelocation = 0;

        s->m_program->enableAttributeArray(vertexlocation);

        s->m_program->setAttributeBuffer(vertexlocation, GL_FLOAT, Vertex2D::positionOffset(), Vertex2D::PositionTupleSize, Vertex2D::stride());

        object.release();
        g->m_vertex.release();
        s->m_program->release();
    }
}


 void GL3DDrawFunctions::DrawModelObject(GL3DWidget * gl, GL3DModel * m, GL3DCamera * camera, QMatrix4x4 ModelMatrix)
{
    for(int i = 0; i < m->GLVAO_List.length(); i++)
    {
        DrawModelGeometryWithMaterialMultipleStart(gl,m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera);
        DrawModelGeometryWithMaterialMultiple(gl,m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera, ModelMatrix);
        DrawModelGeometryWithMaterialMultipleEnd(gl,m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera);
    }
}

 void GL3DDrawFunctions::DrawModelGeometryWithMaterial(GL3DWidget * gl, GL3DGeometry * g, GL3DShader * Shader, QOpenGLVertexArrayObject * vao,GL3DMaterial * mat, GL3DCamera * camera, QMatrix4x4 modelmatrix)
{


}

 void GL3DDrawFunctions::DrawModelGeometryWithMaterialMultipleStart(GL3DWidget * gl, GL3DGeometry * g,GL3DShader * Shader, QOpenGLVertexArrayObject * vao,GL3DMaterial * mat, GL3DCamera * camera)
{
    gl->gl->glDisable(GL_CULL_FACE);

    Shader->m_program->bind();

    Shader->m_program->setUniformValue("Light_Ambient",QVector3D(0.3,0.3,0.3));
    Shader->m_program->setUniformValue("Light_Directional",QVector3D(-1.0,-1.0,0.0));
    Shader->m_program->setUniformValue("Light_Directional_Direction",QVector3D(0.7,0.7,0.7));


    int e = gl->gl->glGetError();

    //set textures
    vao->bind();


    if((mat == 0))
    {
        mat = new GL3DMaterial();
    }

    Shader->m_program->setUniformValue("has_Texture_ka",!(mat->Texture_ka == 0));
    if(!(mat->Texture_ka == 0))
    {
        Shader->ActivateTextureOn(gl,mat->Texture_ka->m_GLTexture,"Texture_ka",0);
    }

    Shader->m_program->setUniformValue("has_Texture_kd",!(mat->Texture_kd == 0));
    if(!(mat->Texture_kd == 0))
    {
        Shader->ActivateTextureOn(gl,mat->Texture_kd->m_GLTexture,"Texture_kd",1);
    }

    Shader->m_program->setUniformValue("has_Texture_ks",!(mat->Texture_ks == 0));
    if(!(mat->Texture_ks == 0))
    {
        Shader->ActivateTextureOn(gl,mat->Texture_ks->m_GLTexture,"Texture_ks",2);
    }

    Shader->m_program->setUniformValue("has_Texture_ns",!(mat->Texture_ns == 0));
    if(!(mat->Texture_ns == 0))
    {
        Shader->ActivateTextureOn(gl,mat->Texture_ns->m_GLTexture,"Texture_ns",3);
    }

    Shader->m_program->setUniformValue("has_Texture_alpha",!(mat->Texture_alpha == 0));
    if(!(mat->Texture_alpha == 0))
    {
        Shader->ActivateTextureOn(gl,mat->Texture_alpha->m_GLTexture,"Texture_alpha",4);
    }

    Shader->m_program->setUniformValue("has_Texture_disp",!(mat->Texture_disp == 0));
    if(!(mat->Texture_disp == 0))
    {
        Shader->ActivateTextureOn(gl,mat->Texture_disp->m_GLTexture,"Texture_disp",5);
    }

    Shader->m_program->setUniformValue("has_Texture_bump",!(mat->Texture_bump == 0));
    if(!(mat->Texture_bump == 0))
    {
        Shader->ActivateTextureOn(gl,mat->Texture_bump->m_GLTexture,"Texture_bump",6);
    }

    Shader->m_program->setUniformValue("has_Texture_normal",!(mat->Texture_normal == 0));
    if(!(mat->Texture_normal == 0))
    {
        Shader->ActivateTextureOn(gl,mat->Texture_normal->m_GLTexture,"Texture_normal",7);
    }

    Shader->m_program->setUniformValue("color_ka",mat->color_ka);
    Shader->m_program->setUniformValue("color_kd",mat->color_kd);
    Shader->m_program->setUniformValue("color_ks",mat->color_ks);
    Shader->m_program->setUniformValue("spec_power",mat->spec_power);
    Shader->m_program->setUniformValue("alpha",(float)mat->alpha);
    Shader->m_program->setUniformValue("illum",mat->illum);

    Shader->m_program->setUniformValue("Light_Ambient",gl->m_World->Light_Ambient);
    Shader->m_program->setUniformValue("Light_Directional",gl->m_World->Light_Directional);
    Shader->m_program->setUniformValue("Light_Direction",gl->m_World->Light_Directional_Direction);

}

 void GL3DDrawFunctions::DrawModelGeometryWithMaterialMultiple(GL3DWidget * gl, GL3DGeometry * g,GL3DShader * Shader, QOpenGLVertexArrayObject * vao,GL3DMaterial * mat, GL3DCamera * camera, QMatrix4x4 ModelMatrix)
{


    Shader->m_program->setUniformValue("CPosition",camera->m_Position);
    Shader->m_program->setUniformValue("Cmatrix",camera->m_CameraMatrix);
    Shader->m_program->setUniformValue("Mmatrix",ModelMatrix);

    gl->gl->glDrawElements(GL_TRIANGLES, g->m_IndexCount, GL_UNSIGNED_INT, 0);

    int e = gl->gl->glGetError();
    if( e != GL_NO_ERROR)
    {
        qDebug() << "opengl error in drawing model object" <<e;
    }


}

 void GL3DDrawFunctions::DrawModelGeometryWithMaterialMultipleEnd(GL3DWidget * gl, GL3DGeometry * g,GL3DShader * Shader, QOpenGLVertexArrayObject * vao,GL3DMaterial * mat, GL3DCamera * camera)
{

    vao->release();

    Shader->m_program->release();

    gl->gl->glEnable(GL_CULL_FACE);
    gl->gl->glCullFace(GL_FRONT);

}

 void GL3DDrawFunctions::DrawModelInstanced(GL3DWidget * gl, GL3DModel * m, GL3DCamera * camera)
{
    for(int i = 0; i < m->GLVAO_List.length(); i++)
    {
        DrawModelGeometryWithMaterialInstanced(gl,m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera, m->m_InstanceCount);

    }

}

 void GL3DDrawFunctions::DrawModelGeometryWithMaterialInstanced(GL3DWidget * gl, GL3DGeometry * g,GL3DShader * Shader, QOpenGLVertexArrayObject * vao,GL3DMaterial * mat, GL3DCamera * camera, int count)
{
    gl->gl->glDisable(GL_CULL_FACE);

    Shader->m_program->bind();

    Shader->m_program->setUniformValue("Light_Ambient",QVector3D(0.3,0.3,0.3));
    Shader->m_program->setUniformValue("Light_Directional",QVector3D(-1.0,-1.0,0.0));
    Shader->m_program->setUniformValue("Light_Directional_Direction",QVector3D(0.7,0.7,0.7));


    //set textures
    vao->bind();


    if((mat == 0))
    {
        mat = new GL3DMaterial();
    }

    Shader->m_program->setUniformValue("has_Texture_ka",!(mat->Texture_ka == 0));
    if(!(mat->Texture_ka == 0))
    {
        Shader->ActivateTextureOn(gl,mat->Texture_ka->m_GLTexture,"Texture_ka",0);
    }

    Shader->m_program->setUniformValue("has_Texture_kd",!(mat->Texture_kd == 0));
    if(!(mat->Texture_kd == 0))
    {
        Shader->ActivateTextureOn(gl,mat->Texture_kd->m_GLTexture,"Texture_kd",1);
    }

    Shader->m_program->setUniformValue("has_Texture_ks",!(mat->Texture_ks == 0));
    if(!(mat->Texture_ks == 0))
    {
        Shader->ActivateTextureOn(gl,mat->Texture_ks->m_GLTexture,"Texture_ks",2);
    }

    Shader->m_program->setUniformValue("has_Texture_ns",!(mat->Texture_ns == 0));
    if(!(mat->Texture_ns == 0))
    {
        Shader->ActivateTextureOn(gl,mat->Texture_ns->m_GLTexture,"Texture_ns",3);
    }

    Shader->m_program->setUniformValue("has_Texture_alpha",!(mat->Texture_alpha == 0));
    if(!(mat->Texture_alpha == 0))
    {
        Shader->ActivateTextureOn(gl,mat->Texture_alpha->m_GLTexture,"Texture_alpha",4);
    }

    Shader->m_program->setUniformValue("has_Texture_disp",!(mat->Texture_disp == 0));
    if(!(mat->Texture_disp == 0))
    {
        Shader->ActivateTextureOn(gl,mat->Texture_disp->m_GLTexture,"Texture_disp",5);
    }

    Shader->m_program->setUniformValue("has_Texture_bump",!(mat->Texture_bump == 0));
    if(!(mat->Texture_bump == 0))
    {
        Shader->ActivateTextureOn(gl,mat->Texture_bump->m_GLTexture,"Texture_bump",6);
    }

    Shader->m_program->setUniformValue("has_Texture_normal",!(mat->Texture_normal == 0));
    if(!(mat->Texture_normal == 0))
    {
        Shader->ActivateTextureOn(gl,mat->Texture_normal->m_GLTexture,"Texture_normal",7);
    }

    Shader->m_program->setUniformValue("color_ka",mat->color_ka);
    Shader->m_program->setUniformValue("color_kd",mat->color_kd);
    Shader->m_program->setUniformValue("color_ks",mat->color_ks);
    Shader->m_program->setUniformValue("spec_power",mat->spec_power);
    Shader->m_program->setUniformValue("alpha",(float)mat->alpha);
    Shader->m_program->setUniformValue("illum",mat->illum);

    Shader->m_program->setUniformValue("Light_Ambient",gl->m_World->Light_Ambient);
    Shader->m_program->setUniformValue("Light_Directional",gl->m_World->Light_Directional);
    Shader->m_program->setUniformValue("Light_Direction",gl->m_World->Light_Directional_Direction);

    Shader->m_program->setUniformValue("CPosition",camera->m_Position);
    Shader->m_program->setUniformValue("Cmatrix",camera->m_CameraMatrix);
    //Shader->m_program->setUniformValue("Mmatrix",ModelMatrix);

    gl->gl->glDrawElementsInstanced(GL_TRIANGLES, g->m_IndexCount, GL_UNSIGNED_INT, 0,count);

    vao->release();

    Shader->m_program->release();

    gl->gl->glEnable(GL_CULL_FACE);
    gl->gl->glCullFace(GL_FRONT);

    int e = gl->gl->glGetError();
    if( e != GL_NO_ERROR)
    {
        qDebug() << "opengl error in drawing model object" <<e;
    }

}


 void GL3DDrawFunctions::DrawModelGLInstanced(GL3DWidget * gl, GL3DModel * m, GL3DCamera * camera, GL3DSurface * surface,double dist_min, double dist_max, double dist_fade, double increment, double rand_location, double rand_rotation, double rand_scale, float * rand)
{
    for(int i = 0; i < m->GLVAO_List.length(); i++)
    {
        DrawModelGeometryWithMaterialGLInstanced(gl,m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera, m->m_InstanceCount, surface,dist_min,dist_max,dist_fade,increment,rand_location,rand_rotation,rand_scale,rand);

    }
}

 void GL3DDrawFunctions::DrawModelGeometryWithMaterialGLInstanced(GL3DWidget * gl, GL3DGeometry * g,GL3DShader * Shader, QOpenGLVertexArrayObject * vao,GL3DMaterial * mat, GL3DCamera * camera, int count, GL3DSurface * surface,double dist_min, double dist_max, double dist_fade, double increment, double rand_location, double rand_rotation, double rand_scale, float * rand)
{
    int rows_display = floor((dist_max + dist_fade)/(increment));
    if((rows_display % 2) == 1)
    {
        rows_display += 1;
    }
    rows_display += 1;
    int cols_display = rows_display;

    int rows_instances = floor(surface->m_ZExtent/increment);
    int cols_instances = floor(surface->m_XExtent/increment);

    qDebug() << increment << rand_location;

    count = rows_display * cols_display;

    Shader->m_program->bind();

    gl->gl->glDisable(GL_CULL_FACE);

    //set textures
    vao->bind();

    Shader->m_program->setUniformValueArray("random",(GLfloat *)rand, 32, 4);

    Shader->m_program->setUniformValue("distances_max_min_fade_increment",QVector4D(dist_max,dist_min,dist_fade,float(increment)));//(float) dist_max);

    Shader->m_program->setUniformValue("dist_min",(float) dist_min);
    Shader->m_program->setUniformValue("dist_max",(float) dist_max);
    Shader->m_program->setUniformValue("dist_fade",(float) dist_fade);
    Shader->m_program->setUniformValue("spacing",(float) increment);

    Shader->m_program->setUniformValue("Cposition",camera->m_Position);

    Shader->m_program->setUniformValue("rand_location",(float) rand_location);
    Shader->m_program->setUniformValue("rand_scale",(float) rand_scale);
    Shader->m_program->setUniformValue("rand_rotation",(float) rand_rotation);
    Shader->m_program->setUniformValue("rows_instances",(float(rows_instances)));
    Shader->m_program->setUniformValue("cols_instances",(float(cols_instances)));
    Shader->m_program->setUniformValue("rows_display",(float(rows_display)));
    Shader->m_program->setUniformValue("cols_display",(float(cols_display)));

    Shader->m_program->setUniformValue("microElevationScaleX",(float) surface->m_CellSize);
    Shader->m_program->setUniformValue("microElevationScaleY",(float) 0.25f);
    Shader->m_program->setUniformValue("microElevationScaleZ",(float) surface->m_CellSize);

    Shader->m_program->setUniformValue("SurfaceExtentX",(float)surface->m_XExtent);
    Shader->m_program->setUniformValue("SurfaceExtentZ",(float)surface->m_ZExtent);

    Shader->m_program->setUniformValue("CellSize",(float)surface->m_CellSize);

    Shader->m_program->setUniformValue("Light_Ambient",gl->m_World->Light_Ambient);
    Shader->m_program->setUniformValue("Light_Directional",gl->m_World->Light_Directional);
    Shader->m_program->setUniformValue("Light_Directional_Direction",gl->m_World->Light_Directional_Direction);

    Shader->ActivateTextureOn(gl,surface->m_Texture,"heightMap",8);
    Shader->ActivateTextureOn(gl,surface->m_Texture_MicroElevation,"microElevation",9);
    Shader->ActivateTextureOn(gl,surface->m_Texture_MicroElevation_Normal,"microElevation_normal",10);
    Shader->ActivateTextureOn(gl,surface->m_Texture_Mask,"Mask",11);

    Shader->ActivateTextureOn(gl,surface->m_Texture_SlopeX,"slopeX",12);
    Shader->ActivateTextureOn(gl,surface->m_Texture_SlopeY,"slopeY",13);

    Shader->ActivateTextureOn(gl,surface->m_Texture_VegCover,"VegCover",14);
    Shader->ActivateTextureOn(gl,surface->m_Texture_VegHeight,"VegHeight",15);
    Shader->ActivateTextureOn(gl,surface->m_Texture_RandomRoughness,"RandomRoughness",16);
    Shader->ActivateTextureOn(gl,surface->m_Texture_Buildings,"Buildings",17);
    Shader->ActivateTextureOn(gl,surface->m_Texture_Roads,"Roads",18);


    if((mat == 0))
    {
        mat = new GL3DMaterial();
    }

    Shader->m_program->setUniformValue("has_Texture_ka",!(mat->Texture_ka == 0));
    if(!(mat->Texture_ka == 0))
    {
        Shader->ActivateTextureOn(gl,mat->Texture_ka->m_GLTexture,"Texture_ka",0);
    }

    Shader->m_program->setUniformValue("has_Texture_kd",!(mat->Texture_kd == 0));
    if(!(mat->Texture_kd == 0))
    {
        Shader->ActivateTextureOn(gl,mat->Texture_kd->m_GLTexture,"Texture_kd",1);
    }

    Shader->m_program->setUniformValue("has_Texture_ks",!(mat->Texture_ks == 0));
    if(!(mat->Texture_ks == 0))
    {
        Shader->ActivateTextureOn(gl,mat->Texture_ks->m_GLTexture,"Texture_ks",2);
    }

    Shader->m_program->setUniformValue("has_Texture_ns",!(mat->Texture_ns == 0));
    if(!(mat->Texture_ns == 0))
    {
        Shader->ActivateTextureOn(gl,mat->Texture_ns->m_GLTexture,"Texture_ns",3);
    }

    Shader->m_program->setUniformValue("has_Texture_alpha",!(mat->Texture_alpha == 0));
    if(!(mat->Texture_alpha == 0))
    {
        Shader->ActivateTextureOn(gl,mat->Texture_alpha->m_GLTexture,"Texture_alpha",4);
    }

    Shader->m_program->setUniformValue("has_Texture_disp",!(mat->Texture_disp == 0));
    if(!(mat->Texture_disp == 0))
    {
        Shader->ActivateTextureOn(gl,mat->Texture_disp->m_GLTexture,"Texture_disp",5);
    }

    Shader->m_program->setUniformValue("has_Texture_bump",!(mat->Texture_bump == 0));
    if(!(mat->Texture_bump == 0))
    {
        Shader->ActivateTextureOn(gl,mat->Texture_bump->m_GLTexture,"Texture_bump",6);
    }

    Shader->m_program->setUniformValue("has_Texture_normal",!(mat->Texture_normal == 0));
    if(!(mat->Texture_normal == 0))
    {
        Shader->ActivateTextureOn(gl,mat->Texture_normal->m_GLTexture,"Texture_normal",7);
    }

    Shader->m_program->setUniformValue("color_ka",mat->color_ka);
    Shader->m_program->setUniformValue("color_kd",mat->color_kd);
    Shader->m_program->setUniformValue("color_ks",mat->color_ks);
    Shader->m_program->setUniformValue("spec_power",mat->spec_power);
    Shader->m_program->setUniformValue("alpha",(float)mat->alpha);
    Shader->m_program->setUniformValue("illum",mat->illum);

    Shader->m_program->setUniformValue("Cmatrix",camera->m_CameraMatrix);
    //Shader->m_program->setUniformValue("Mmatrix",ModelMatrix);

    gl->gl->glDrawElementsInstanced(GL_TRIANGLES, g->m_IndexCount, GL_UNSIGNED_INT, 0,count);

    vao->release();

    Shader->m_program->release();

    gl->gl->glEnable(GL_CULL_FACE);
    gl->gl->glCullFace(GL_FRONT);

    int e = gl->gl->glGetError();
    if( e != GL_NO_ERROR)
    {
        qDebug() << "opengl error in drawing model object" <<e;
    }

}
