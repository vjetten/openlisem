/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Author: Victor Jetten
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/

#include <3D/Graphics/GL3DTextures.h>

void GL3DTextures::Create(GL3DWidget * widget)
{

    m_Widget = widget;
    m_LoadedTextureList.clear();

}

GL3DTexture * GL3DTextures::LoadTextureFromMap(bool debug, cTMap * elevation, int res_x,int res_y, bool data,bool mask,bool fill, GL3DColorRamp * color_ramp)
{
    GL3DTexture * t = new GL3DTexture();

    t->CreateTexture(debug, m_Widget,elevation,res_x,res_y,data,mask,color_ramp);

    this->m_LoadedTextureList.append(t);

    return t;
}

GL3DTexture * GL3DTextures::LoadTextureFromPerlin(int lx, int ly,QList<float> wave, QList<float> amp)
{
    GL3DTexture * t = new GL3DTexture();

    t->CreateTexturePerlin(m_Widget,lx,ly,wave,amp);

    this->m_LoadedTextureList.append(t);

    return t;
}

GL3DTexture * GL3DTextures::LoadTextureFromMatrix(float * ,int lx, int ly)
{

    return 0;
}

GL3DTexture * GL3DTextures::LoadTextureFromFile(QString file,bool add, bool is_repeat)
{
    GL3DTexture * t = new GL3DTexture();

    t->CreateTexture(m_Widget,file,is_repeat);
    if(add)
    {
        this->m_LoadedTextureList.append(t);
    }
    return t;

}

GL3DTexture * GL3DTextures::LoadCubeTextureFromFile(QString file,bool add)
{
    GL3DTexture * t = new GL3DTexture();

    t->CreateCubeTexture(m_Widget,file);
    if(add)
    {
        this->m_LoadedTextureList.append(t);
    }
    return t;
}

void GL3DTextures::ClearUnused()
{


}

void GL3DTextures::Destroy()
{



}


void GL3DTextures::CreateFrameBuffers(GL3DWidget * widget,int w, int h)
{
    qDebug() << "create frame buffers " << w << h << "gl context " << widget << widget->gl;

    if(w == 0 || h == 0)
    {
        return;
    }
    if(this->bufferscreated)
    {
        widget->gl->glDeleteTextures(1,&RenderTexture);
        widget->gl->glDeleteTextures(1,&RenderTextureCopy);
        widget->gl->glDeleteTextures(1,&RenderTextureWater);
        widget->gl->glDeleteTextures(1,&LocationTexture);
        widget->gl->glDeleteTextures(1,&NormalTexture);
        widget->gl->glDeleteTextures(1,&InfoTexture);
        widget->gl->glDeleteTextures(1,&ChannelTexture);
        widget->gl->glDeleteTextures(1,&ChannelInfoTexture);

        //weird, removing frame buffers gives error when not in debugging mode
        /*qDebug() << "deleted" << FramebufferCopy << ChannelFramebuffer << Framebuffer << FramebufferWater;
        widget->gl->glDeleteFramebuffers(GL_FRAMEBUFFER, &FramebufferCopy);
        qDebug() << "deleted";
        widget->gl->glDeleteFramebuffers(GL_FRAMEBUFFER, &ChannelFramebuffer);
        qDebug() << "deleted";
        widget->gl->glDeleteFramebuffers(GL_FRAMEBUFFER, &Framebuffer);
        qDebug() << "deleted";
        widget->gl->glDeleteFramebuffers(GL_FRAMEBUFFER, &FramebufferWater);*/

    }

    widget->gl->glGenFramebuffers(1, &Framebuffer);
    widget->gl->glBindFramebuffer(GL_FRAMEBUFFER, Framebuffer);

    widget->gl->glGenTextures(1, &RenderTexture);
    widget->gl->glGenTextures(1, &RenderTextureCopy);
    widget->gl->glGenTextures(1, &LocationTexture);
    widget->gl->glGenTextures(1, &NormalTexture);
    widget->gl->glGenTextures(1, &InfoTexture);

    widget->gl->glBindTexture(GL_TEXTURE_2D, RenderTexture);
    widget->gl->glTexImage2D(GL_TEXTURE_2D, 0,GL_RGBA32F, w, h, 0,GL_RGBA, GL_FLOAT, 0);
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    widget->gl->glBindTexture(GL_TEXTURE_2D, RenderTextureCopy);
    widget->gl->glTexImage2D(GL_TEXTURE_2D, 0,GL_RGBA32F, w, h, 0,GL_RGBA, GL_FLOAT, 0);
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    widget->gl->glBindTexture(GL_TEXTURE_2D, LocationTexture);
    widget->gl->glTexImage2D(GL_TEXTURE_2D, 0,GL_RGBA32F, w, h, 0,GL_RGBA, GL_FLOAT, 0);
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    widget->gl->glBindTexture(GL_TEXTURE_2D, NormalTexture);
    widget->gl->glTexImage2D(GL_TEXTURE_2D, 0,GL_RGBA32F, w, h, 0,GL_RGBA, GL_FLOAT, 0);
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    widget->gl->glBindTexture(GL_TEXTURE_2D, InfoTexture);
    widget->gl->glTexImage2D(GL_TEXTURE_2D, 0,GL_RGBA32F, w, h, 0,GL_RGBA, GL_FLOAT, 0);
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);


    widget->gl->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, RenderTexture, 0);
    widget->gl->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, LocationTexture, 0);
    widget->gl->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, NormalTexture, 0);
    widget->gl->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT3, InfoTexture, 0);

    widget->gl->glGenRenderbuffers(1, &depthrenderbuffer);
    widget->gl->glBindRenderbuffer(GL_RENDERBUFFER, depthrenderbuffer);
    widget->gl->glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, w,h);
    widget->gl->glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthrenderbuffer);


    widget->gl->glGenFramebuffers(1, &FramebufferCopy);
    widget->gl->glBindFramebuffer(GL_FRAMEBUFFER, FramebufferCopy);
    widget->gl->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT4, RenderTextureCopy, 0);

    widget->gl->glGenTextures(1, &RenderTextureCopy);

    widget->gl->glBindTexture(GL_TEXTURE_2D, RenderTextureCopy);
    widget->gl->glTexImage2D(GL_TEXTURE_2D, 0,GL_RGBA32F, w, h, 0,GL_RGBA, GL_FLOAT, 0);
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    widget->gl->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, RenderTextureCopy, 0);


    widget->gl->glGenFramebuffers(1, &ChannelFramebuffer);
    widget->gl->glBindFramebuffer(GL_FRAMEBUFFER, ChannelFramebuffer);

    widget->gl->glGenTextures(1, &ChannelTexture);
    widget->gl->glGenTextures(1, &ChannelInfoTexture);

    widget->gl->glBindTexture(GL_TEXTURE_2D, ChannelTexture);
    widget->gl->glTexImage2D(GL_TEXTURE_2D, 0,GL_RGBA32F, w, h, 0,GL_RGBA, GL_FLOAT, 0);
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    widget->gl->glBindTexture(GL_TEXTURE_2D, ChannelInfoTexture);
    widget->gl->glTexImage2D(GL_TEXTURE_2D, 0,GL_RGBA32F, w, h, 0,GL_RGBA, GL_FLOAT, 0);
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    widget->gl->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, ChannelTexture, 0);
    widget->gl->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, ChannelInfoTexture, 0);

    widget->gl->glGenRenderbuffers(1, &channeldepthrenderbuffer);
    widget->gl->glBindRenderbuffer(GL_RENDERBUFFER, channeldepthrenderbuffer);
    widget->gl->glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, w,h);
    widget->gl->glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, channeldepthrenderbuffer);



    widget->gl->glGenFramebuffers(1, &FramebufferWater);
    widget->gl->glBindFramebuffer(GL_FRAMEBUFFER, FramebufferWater);

    widget->gl->glGenTextures(1, &RenderTextureWater);

    widget->gl->glBindTexture(GL_TEXTURE_2D, RenderTextureWater);
    widget->gl->glTexImage2D(GL_TEXTURE_2D, 0,GL_RGBA32F, w, h, 0,GL_RGBA, GL_FLOAT, 0);
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    widget->gl->glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, RenderTextureWater, 0);

    widget->gl->glGenRenderbuffers(1, &depthrenderbuffer);
    widget->gl->glBindRenderbuffer(GL_RENDERBUFFER, depthrenderbuffer);
    widget->gl->glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, w,h);
    widget->gl->glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthrenderbuffer);


    bufferscreated = true;

}

void GL3DTexture::CreateTexture(GL3DWidget * widget,float * ,int lx, int ly)
{




}

void GL3DTexture::CreateTexturePerlin(GL3DWidget *widget, int lx, int ly, QList<float> wave, QList<float> amp)
{

    int n = lx * ly;
    GLfloat*data = (float *)malloc(n * sizeof(GLfloat) );
    for(int y = 0; y < ly;y++)
    {
        for(int x = 0; x < lx;x++)
        {
            data[y * lx + x] = 0;
        }
    }
    for(int i = 0; i < wave.length();i++)
    {

        PerlinNoise *p = new PerlinNoise(i*123456);
        for(int y = 0; y < ly;y++)
        {
            for(int x = 0; x < lx;x++)
            {
                double rx = float(x) * 1.0/wave.at(i);
                double ry = float(y) * 1.0/wave.at(i);
                double val = p->noise(rx,ry,0);

                data[y * lx + x] += val * amp.at(i);

            }

        }
        delete p;
    }

    //creat a texture and store data
    widget->gl->glGenTextures(1, &m_GLTexture);
    widget->gl->glActiveTexture(GL_TEXTURE0);
    widget->gl->glBindTexture(GL_TEXTURE_2D,  m_GLTexture);
    widget->gl->glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, lx,ly, 0, GL_RED, GL_FLOAT, data);

    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    widget->gl->glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
    widget->gl->glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    free(data);

    int e = widget->gl->glGetError();
    if( e != GL_NO_ERROR)
    {
        qDebug() << "opengl error in creating perlin noise texture" << e;
    }

    this->is_created = true;

}

void GL3DTexture::CreateTexture(bool debug, GL3DWidget * widget,cTMap * elevation, int res_x,int res_y, bool data,bool mask,bool fill, GL3DColorRamp * color_ramp )
{
    double emax = 0;
    double emin = 0;
    bool first = true;

    FOR_ROW_COL_MV(elevation,elevation)
    {
        if(first)
        {
            first = false;
            emax = elevation->Drc;
            emin = elevation->Drc;
        }else
        {

            emin = std::min(emin,elevation->Drc);
            emax = std::max(emax,elevation->Drc);
        }
    }

    if(!data)
    {
        this->m_QImage = QImage(elevation->nrCols(), elevation->nrRows(), QImage::Format_ARGB32);

        QVector4D value;

        for(int x = 0; x < elevation->nrCols();x++)
        {
            for(int y = 0; y < elevation->nrRows();y++)
            {
                if(pcr::isMV(elevation->data[y][x]))
                {
                    value = color_ramp->GetColor(emin,emax,elevation->data[y][x]);
                    m_QImage.setPixel(x,y, (QColor(value.x(),value.y(),value.z(),value.w()).rgba()));
                }else
                {
                    m_QImage.setPixel(x,y, (QColor(0,0,0,0).rgba()));
                }
            }
        }


        m_QGLImage = QGLWidget::convertToGLFormat(m_QImage);

        widget->gl->glGenTextures(1,&m_GLTexture);
        widget->gl->glActiveTexture(GL_TEXTURE0);
        widget->gl->glBindTexture( GL_TEXTURE_2D, m_GLTexture );
        widget->gl->glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, m_QGLImage.width(), m_QGLImage.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, m_QGLImage.bits());

        widget->gl->glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        widget->gl->glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        widget->gl->glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
        widget->gl->glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
        widget->gl->glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );

    }else if(!mask)
    {
        int n = elevation->nrCols()*elevation->nrRows();
        GLfloat*data = (float *)malloc(n * sizeof(GLfloat) );

        int npos = 0;

        for(int x = 0; x < elevation->nrCols();x++)
        {
            for(int y = 0; y < elevation->nrRows();y++)
            {
                float elevation_c = 0;
                int count = 0;
                if(!pcr::isMV(elevation->data[y][x]))
                {
                    if(debug)
                    {
                        elevation_c = elevation->data[y][x];
                        qDebug() <<x << y << elevation_c;
                    }else
                    {
                        elevation_c = elevation->data[y][x];
                    }
                }else
                {
                    elevation_c = emin;
                }

                data[(y*elevation->nrCols() + x)]=elevation_c;
            }
        }

        //creat a texture and store data
        widget->gl->glGenTextures(1, &m_GLTexture);
        widget->gl->glActiveTexture(GL_TEXTURE0);
        widget->gl->glBindTexture(GL_TEXTURE_2D,  m_GLTexture);
        widget->gl->glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, elevation->nrCols(),elevation->nrRows(), 0, GL_RED, GL_FLOAT, data);

        widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
        widget->gl->glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
        widget->gl->glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        free(data);
    }else
    {
        int n = elevation->nrCols()*elevation->nrRows();
        GLfloat*data = (GLfloat *)malloc(n * sizeof(GLfloat) );

        for(int x = 0; x < elevation->nrCols();x++)
        {
            for(int y = 0; y < elevation->nrRows();y++)
            {
                float elevation_c = 0;
                int count = 0;
                if(!pcr::isMV(elevation->data[y][x]))
                {
                    data[(y*elevation->nrCols() + x)]= (GLfloat) ( 1.0);

                }else
                {
                    data[(y*elevation->nrCols() + x)]= (GLfloat) 0.0;
                }
            }
        }

        //creat a texture and store data
        widget->gl->glGenTextures(1, &m_GLTexture);
        widget->gl->glActiveTexture(GL_TEXTURE0);
        widget->gl->glBindTexture(GL_TEXTURE_2D,  m_GLTexture);
        widget->gl->glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, elevation->nrCols(),elevation->nrRows(), 0, GL_RED, GL_FLOAT, data);

        widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
        widget->gl->glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
        widget->gl->glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

        free(data);
    }

    int e = widget->gl->glGetError();
    if( e != GL_NO_ERROR)
    {
        qDebug() << "opengl error in converting map to texture" << e;
    }

    this->is_created = true;
}

void GL3DTexture::CreateTextureDirectPath(GL3DWidget * widget,QString file, bool is_repeat)
{
    QString long_name = file;
    bool load = m_QImage.load(long_name);
    if(!load)
    {
        qDebug() << "resource not available" << long_name;
        return;
    }

    this->width = m_QImage.width();
    this->height = m_QImage.height();
    m_QGLImage = QGLWidget::convertToGLFormat(m_QImage);

    widget->gl->glGenTextures(1,&m_GLTexture);
    widget->gl->glActiveTexture(GL_TEXTURE0);
    widget->gl->glBindTexture( GL_TEXTURE_2D, m_GLTexture );
    widget->gl->glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, m_QGLImage.width(), m_QGLImage.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, m_QGLImage.bits());
    widget->gl->glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
    widget->gl->glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
    if(is_repeat)
    {
        widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    }else
    {
        widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    }
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
    widget->gl->glGenerateMipmap(GL_TEXTURE_2D);
    this->is_created = true;

}

void GL3DTexture::CreateTexture(GL3DWidget * widget,QString file, bool is_repeat)
{
    QString long_name = widget->m_Directory + "/" + GL3D_DIR_TEXTURES + file;
    CreateTextureDirectPath(widget,long_name,is_repeat);


}
void GL3DTexture::CreateCubeTexture(GL3DWidget * widget,QString file)
{
    QString long_name = widget->m_Directory + "/" +GL3D_DIR_TEXTURES + file;
    bool load = m_QImage.load(long_name);
    if(!load)
    {
        qDebug() << "resource not available" << long_name;
        return;
    }

    this->width = m_QImage.width();
    this->height = m_QImage.height();

    int width_face = this->width/4;
    int height_face = this->height/3;

    QRect rect_pos_x = QRect(2 * width_face+1,1 * height_face+1,width_face-2,height_face-2);
    QRect rect_neg_x = QRect(0 * width_face+1,1 * height_face+1,width_face-2,height_face-2);
    QRect rect_pos_y = QRect(1 * width_face+1,0 * height_face+1,width_face-2,height_face-2);
    QRect rect_neg_y = QRect(1 * width_face+1,2 * height_face+1,width_face-2,height_face-2);
    QRect rect_neg_z = QRect(3 * width_face+1,1 * height_face+1,width_face-2,height_face-2);
    QRect rect_pos_z = QRect(1 * width_face+1,1 * height_face+1,width_face-2,height_face-2);

    i_pos_x = createSubImage(&m_QImage,rect_pos_x).mirrored();
    i_neg_x = createSubImage(&m_QImage,rect_neg_x).mirrored();
    i_pos_y = createSubImage(&m_QImage,rect_pos_y).mirrored();
    i_neg_y = createSubImage(&m_QImage,rect_neg_y).mirrored();
    i_pos_z = createSubImage(&m_QImage,rect_pos_z).mirrored();
    i_neg_z = createSubImage(&m_QImage,rect_neg_z).mirrored();

    m_QGLImage_pos_x = QGLWidget::convertToGLFormat(i_pos_x);
    m_QGLImage_neg_x = QGLWidget::convertToGLFormat(i_neg_x);
    m_QGLImage_pos_y = QGLWidget::convertToGLFormat(i_pos_y);
    m_QGLImage_neg_y = QGLWidget::convertToGLFormat(i_neg_y);
    m_QGLImage_pos_z = QGLWidget::convertToGLFormat(i_pos_z);
    m_QGLImage_neg_z = QGLWidget::convertToGLFormat(i_neg_z);

    widget->gl->glEnable(GL_TEXTURE_CUBE_MAP);
    widget->gl->glGenTextures(1,&m_GLTexture);
    widget->gl->glActiveTexture(GL_TEXTURE0);
    widget->gl->glBindTexture( GL_TEXTURE_CUBE_MAP, m_GLTexture );

    widget->gl->glTexImage2D( GL_TEXTURE_CUBE_MAP_POSITIVE_X, 0, GL_RGBA, m_QGLImage_pos_x.width(), m_QGLImage_pos_x.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, m_QGLImage_pos_x.bits());
    widget->gl->glTexImage2D( GL_TEXTURE_CUBE_MAP_NEGATIVE_X, 0, GL_RGBA, m_QGLImage_neg_x.width(), m_QGLImage_neg_x.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, m_QGLImage_neg_x.bits());
    widget->gl->glTexImage2D( GL_TEXTURE_CUBE_MAP_POSITIVE_Y, 0, GL_RGBA, m_QGLImage_pos_y.width(), m_QGLImage_pos_y.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, m_QGLImage_pos_y.bits());
    widget->gl->glTexImage2D( GL_TEXTURE_CUBE_MAP_NEGATIVE_Y, 0, GL_RGBA, m_QGLImage_neg_y.width(), m_QGLImage_neg_y.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, m_QGLImage_neg_y.bits());
    widget->gl->glTexImage2D( GL_TEXTURE_CUBE_MAP_POSITIVE_Z, 0, GL_RGBA, m_QGLImage_pos_z.width(), m_QGLImage_pos_z.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, m_QGLImage_pos_z.bits());
    widget->gl->glTexImage2D( GL_TEXTURE_CUBE_MAP_NEGATIVE_Z, 0, GL_RGBA, m_QGLImage_neg_z.width(), m_QGLImage_neg_z.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, m_QGLImage_neg_z.bits());

    widget->gl->glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    widget->gl->glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    widget->gl->glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    widget->gl->glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    widget->gl->glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

    widget->gl->glBindTexture( GL_TEXTURE_CUBE_MAP, 0 );

    this->is_created = true;
    widget->gl->glDisable(GL_TEXTURE_CUBE_MAP);
}

//needs original image, just point to right data
QImage GL3DTexture::createSubImage(QImage* image, const QRect & rect) {
    size_t offset = rect.x() * image->depth() / 8
                    + rect.y() * image->bytesPerLine();
    return QImage(image->bits() + offset, rect.width(), rect.height(),
                  image->bytesPerLine(), image->format());
}

void GL3DTexture::UpdateTextureFromMap(GL3DWidget * widget, cTMap * elevation, int res_x,int res_y, bool data, bool mask, bool fill, GL3DColorRamp * color_ramp)
{

    double emax = 0;
    double emin = 0;
    bool first = true;

    FOR_ROW_COL_MV(elevation,elevation)
    {
        if(first)
        {
            first = false;
            emax = elevation->Drc;
            emin = elevation->Drc;
        }else
        {

            emin = std::min(emin,elevation->Drc);
            emax = std::max(emax,elevation->Drc);
        }
    }

    if(!data)
    {
        this->m_QImage = QImage(elevation->nrCols(), elevation->nrRows(), QImage::Format_ARGB32);

        QVector4D value;

        for(int x = 0; x < elevation->nrCols();x++)
        {
            for(int y = 0; y < elevation->nrRows();y++)
            {
                if(pcr::isMV(elevation->data[y][x]))
                {
                    value = color_ramp->GetColor(emin,emax,elevation->data[y][x]);
                    m_QImage.setPixel(x,y, (QColor(value.x(),value.y(),value.z(),value.w()).rgba()));
                }else
                {
                    m_QImage.setPixel(x,y, (QColor(0,0,0,0).rgba()));
                }
            }
        }


        m_QGLImage = QGLWidget::convertToGLFormat(m_QImage);

        widget->gl->glActiveTexture(GL_TEXTURE0);
        widget->gl->glBindTexture( GL_TEXTURE_2D, m_GLTexture );
        widget->gl->glTexSubImage2D( GL_TEXTURE_2D, 0,0,0 , m_QGLImage.width(), m_QGLImage.height(), GL_RGBA, GL_UNSIGNED_BYTE, m_QGLImage.bits());

    }else if(!mask)
    {
        int n = elevation->nrCols()*elevation->nrRows();
        GLfloat*data = (float *)malloc(n * sizeof(GLfloat) );

        int npos = 0;

        for(int x = 0; x < elevation->nrCols();x++)
        {
            for(int y = 0; y < elevation->nrRows();y++)
            {
                float elevation_c = 0;
                int count = 0;
                if(!pcr::isMV(elevation->data[y][x]))
                {
                    elevation_c = elevation->data[y][x];
                }else
                {
                    elevation_c = emin;
                }

                data[(y*elevation->nrCols() + x)]=elevation_c;
            }
        }

        //creat a texture and store data
        widget->gl->glActiveTexture(GL_TEXTURE0);
        widget->gl->glBindTexture(GL_TEXTURE_2D,  m_GLTexture);
        widget->gl->glTexSubImage2D(GL_TEXTURE_2D, 0,0, 0, elevation->nrCols(),elevation->nrRows(), GL_RED, GL_FLOAT, data);

        free(data);
    }else
    {
        int n = elevation->nrCols()*elevation->nrRows();
        GLfloat*data = (GLfloat *)malloc(n * sizeof(GLfloat) );

        for(int x = 0; x < elevation->nrCols();x++)
        {
            for(int y = 0; y < elevation->nrRows();y++)
            {
                float elevation_c = 0;
                int count = 0;
                if(!pcr::isMV(elevation->data[y][x]))
                {
                    data[(y*elevation->nrCols() + x)]= (GLfloat) ( 1.0);

                }else
                {
                    data[(y*elevation->nrCols() + x)]= (GLfloat) 0.0;
                }
            }
        }

        //creat a texture and store data
        widget->gl->glActiveTexture(GL_TEXTURE0);
        widget->gl->glBindTexture(GL_TEXTURE_2D,  m_GLTexture);
        widget->gl->glTexSubImage2D(GL_TEXTURE_2D, 0,0,0, elevation->nrCols(),elevation->nrRows(), GL_RED, GL_FLOAT, data);

        free(data);
    }

    int e = widget->gl->glGetError();
    if( e != GL_NO_ERROR)
    {
        qDebug() << "opengl error in converting map to texture" << e;
    }

    this->is_created = true;



}


void GL3DTexture::ClearTexture(GL3DWidget * widget)
{
    widget->gl->glDeleteTextures(1,&(this->m_GLTexture));
    is_created = false;


}
