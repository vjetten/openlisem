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

#include <3D/Graphics/GL3DModels.h>
#include <QStringList>

void GL3DModels::Create(GL3DWidget * widget)
{
    m_Widget= widget;


}

void GL3DModels::ClearUnused()
{


}

void GL3DModels::Destroy()
{


}

QList<GL3DMaterial*> GL3DModel::LoadMaterialsFile(GL3DWidget * widget,QString in_file, QString path)
{
    QList<GL3DMaterial *> temp_list;

    GL3DMaterial * current_material = 0;

    //input is long filename
    QFile file(in_file);

    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return temp_list ;

    QTextStream in(&file);
    while (!in.atEnd())
    {

        QString line = in.readLine();
        line = line.trimmed();

        QStringList list = line.split(QRegExp("\\s+"),QString::SplitBehavior::SkipEmptyParts );

        QStringList options;

        for(int i = 0; i < list.length(); i++)
        {
            if(list.at(i).trimmed().startsWith(QChar('-')))
            {
                //option detected, not supported
                options.append(list.at(i));

                //remove from list and decrease i, since indexes of list have changed
                list.removeAt(i); i--;

            }
        }
        line = QString("");
        for(int i = 0; i < list.length(); i++)
        {
            line = line + " " +  list.at(i);
        }

        //empty line
        if(list.length() == 0)
        {
            continue;
        }

        QString first = list.at(0);

        //comment
        if(first.at(0) == QLatin1Char('#'))
        {
            continue;
        }

        if(first.compare(QString("newmtl")) == 0)
        {
            current_material = new GL3DMaterial();
            temp_list.append(current_material);
            current_material->name = line.right(line.length() - (first.length() + 1));

        }else if(first.compare(QString("Ka")) == 0)
        {
            vec3d temp_vec;
            sscanf(line.right(line.length() - (first.length() + 1)).toLatin1(), "%f %f %f", &temp_vec.x, &temp_vec.y, &temp_vec.z );
            if(current_material != 0)
            {
                current_material->color_ka = QVector3D(temp_vec.x,temp_vec.y,temp_vec.z);
            }

        }else if(first.compare(QString("Kd")) == 0)
        {
            vec3d temp_vec;
            sscanf(line.right(line.length() - (first.length() + 1)).toLatin1(), "%f %f %f", &temp_vec.x, &temp_vec.y, &temp_vec.z );
            if(current_material != 0)
            {
                current_material->color_kd = QVector3D(temp_vec.x,temp_vec.y,temp_vec.z);
            }

        }else if(first.compare(QString("Ks")) == 0)
        {
            vec3d temp_vec;
            sscanf(line.right(line.length() - (first.length() + 1)).toLatin1(), "%f %f %f", &temp_vec.x, &temp_vec.y, &temp_vec.z );
            if(current_material != 0)
            {
                current_material->color_ks = QVector3D(temp_vec.x,temp_vec.y,temp_vec.z);
            }

        }else if(first.compare(QString("Ns")) == 0)
        {
            vec3d temp_vec;
            sscanf(line.right(line.length() - (first.length() + 1)).toLatin1(), "%f %f %f", &temp_vec.x, &temp_vec.y, &temp_vec.z );
            if(current_material != 0)
            {
                current_material->color_ks = QVector3D(temp_vec.x,temp_vec.y,temp_vec.z);
            }

        }
        if(first.compare(QString("d")) == 0)
        {
            vec3d temp_vec;
            sscanf(line.right(line.length() - (first.length() + 1)).toLatin1(), "%f", &temp_vec.x);
            if(current_material != 0)
            {
                current_material->alpha = (temp_vec.x);
            }

        }else if(first.compare(QString("Tr")) == 0)
        {
            vec3d temp_vec;
            sscanf(line.right(line.length() - (first.length() + 1)).toLatin1(), "%f", &temp_vec.x);
            if(current_material != 0)
            {
                current_material->alpha =  1.0 -(temp_vec.x);
            }
        }
        if(first.compare(QString("map_Ka")) == 0)
        {
            if(current_material != 0)
            {
                current_material->Texture_ka = new GL3DTexture();
                current_material->Texture_ka->CreateTextureDirectPath(widget,path + line.right(line.length() - (first.length() + 1)),true);
            }
        }else if(first.compare(QString("map_Kd")) == 0)
        {
            if(current_material != 0)
            {
                current_material->Texture_kd = new GL3DTexture();
                current_material->Texture_kd->CreateTextureDirectPath(widget,path + line.right(line.length() - (first.length() + 1)),true);
            }
        }else if(first.compare(QString("map_Ks")) == 0)
        {
            if(current_material != 0)
            {
                current_material->Texture_ks = new GL3DTexture();
                current_material->Texture_ks->CreateTextureDirectPath(widget,path + line.right(line.length() - (first.length() + 1)),true);
            }
        }else if(first.compare(QString("map_Ns")) == 0)
        {
            if(current_material != 0)
            {
                current_material->Texture_ns = new GL3DTexture();
                current_material->Texture_ns->CreateTextureDirectPath(widget,path + line.right(line.length() - (first.length() + 1)),true);
            }

        }else if(first.compare(QString("map_d")) == 0)
        {
            if(current_material != 0)
            {
                current_material->Texture_alpha = new GL3DTexture();
                current_material->Texture_alpha->CreateTextureDirectPath(widget,path + line.right(line.length() - (first.length() + 1)),true);
            }

        }else if(first.compare(QString("map_bump")) == 0 || first.compare(QString("bump")) == 0)
        {
            if(current_material != 0)
            {
                current_material->Texture_bump = new GL3DTexture();
                current_material->Texture_bump->CreateTextureDirectPath(widget,path + line.right(line.length() - (first.length() + 1)),true);
            }
        }else if(first.compare(QString("map_disp")) == 0 || first.compare(QString("disp")) == 0)
        {
            if(current_material != 0)
            {
                current_material->Texture_disp = new GL3DTexture();
                current_material->Texture_disp->CreateTextureDirectPath(widget,path + line.right(line.length() - (first.length() + 1)),true);
            }
        }else if(first.compare(QString("map_normal")) == 0 || first.compare(QString("normal")) == 0)
        {
            if(current_material != 0)
            {
                current_material->Texture_normal = new GL3DTexture();
                current_material->Texture_normal->CreateTextureDirectPath(widget,path + line.right(line.length() - (first.length() + 1)),true);
            }
        }else if(first.compare(QString("illum")) == 0 )
        {
            if(current_material != 0)
            {
                current_material->Texture_normal = new GL3DTexture();
                current_material->Texture_normal->CreateTextureDirectPath(widget,path + line.right(line.length() - (first.length() + 1)),true);
            }
        }

    }

    return temp_list;
}



void GL3DModel::LoadObjectFile(GL3DWidget * widget,QString in_file)
{

    Destroy(widget);

    std::vector< unsigned int > vertexIndices, uvIndices, normalIndices;
    std::vector< QVector3D > temp_vertices;
    std::vector< QVector2D > temp_uvs;
    std::vector< QVector3D > temp_normals;

    bool finish_geometry = false;
    bool unfinished_data = false;
    bool materials_loaded = false;
    bool set_new_material = false;
    QString new_material_name;
    int current_material = 0;

    GL3DMaterial *m_null = new GL3DMaterial();

    Material_List.append(m_null);

    QFile file(widget->m_Directory + "/" + GL3D_DIR_SHADERS + in_file);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return;

    QString path = QFileInfo(widget->m_Directory + "/" + GL3D_DIR_SHADERS + in_file).absolutePath();

    QTextStream in(&file);
    while (!in.atEnd()) {

        QString line = in.readLine();
        line = line.trimmed();

        QStringList list = line.split(QRegExp("\\s+"),QString::SplitBehavior::SkipEmptyParts );

        //empty line
        if(list.length() == 0)
        {
            continue;
        }

        QString first = list.at(0);

        //comment
        if(first.at(0) == QLatin1Char('#'))
        {
            continue;
        }

        //read the material library
        if(first.compare(QString("mtllib")) == 0)
        {

            Material_List.append(LoadMaterialsFile(widget,path + "/" + line.right(line.length() - (first.length() + 1)),path));

            materials_loaded = true;

        }else if(first.compare(QString("v")) == 0)
        {
            vec3d temp_vec;
            sscanf(line.right(line.length() - (first.length() + 1)).toLatin1(), "%f %f %f", &temp_vec.x, &temp_vec.y, &temp_vec.z );
            temp_vertices.push_back(QVector3D(temp_vec.x,temp_vec.y,temp_vec.z));

        }else if(first.compare(QString("vt")) == 0)
        {
            vec3d temp_vec;
            sscanf(line.right(line.length() - (first.length() + 1)).toLatin1(), "%f %f", &temp_vec.x, &temp_vec.y );
            temp_uvs.push_back(QVector2D(temp_vec.x,temp_vec.y));

        }else if(first.compare(QString("vn")) == 0)
        {
            vec3d temp_vec;
            sscanf(line.right(line.length() - (first.length() + 1)).toLatin1(), "%f %f %f", &temp_vec.x, &temp_vec.y, &temp_vec.z );
            temp_normals.push_back(QVector3D(temp_vec.x,temp_vec.y,temp_vec.z));

        //new object
        }else if(first.compare(QString("o")) == 0)
        {

            finish_geometry = true;

        }else if(first.compare(QString("g")) == 0)
        {

            finish_geometry = true;


        //use a material for the following geometry
        }else if(first.compare(QString("usemtl")) == 0)
        {

            finish_geometry = true;
            set_new_material = true;
            new_material_name = line.right(line.length() - first.length() + 1);


        }else if(first.compare(QString("f")) == 0)
        {
            unsigned int vertexIndex[3], uvIndex[3], normalIndex[3];

            QByteArray oString = line.right(line.length() - first.length() + 1).toLatin1();
            const char *pszString = oString.constData();

            int matches = sscanf(pszString, "%d/%d/%d %d/%d/%d %d/%d/%d", &vertexIndex[0], &uvIndex[0], &normalIndex[0], &vertexIndex[1], &uvIndex[1], &normalIndex[1], &vertexIndex[2], &uvIndex[2], &normalIndex[2] );
            if (matches != 9){
                printf(".obj model file must only contain triangles");
                return;
            }

            vertexIndices.push_back(vertexIndex[0]);
            vertexIndices.push_back(vertexIndex[1]);
            vertexIndices.push_back(vertexIndex[2]);
            uvIndices.push_back(uvIndex[0]);
            uvIndices.push_back(uvIndex[1]);
            uvIndices.push_back(uvIndex[2]);
            normalIndices.push_back(normalIndex[0]);
            normalIndices.push_back(normalIndex[1]);
            normalIndices.push_back(normalIndex[2]);

            unfinished_data = true;

        }

        if(finish_geometry)
        {
            if(unfinished_data)
            {
                GL3DGeometry * g = new GL3DGeometry();
                int nr_vertex = vertexIndices.size()/3;

                Vertex * vertexdata = (Vertex*) malloc(sizeof(Vertex) * nr_vertex * 3);
                GLuint * indexdata =(GLuint*)  malloc(sizeof(int) * nr_vertex * 3);

                for(int i = 0; i < nr_vertex; i++)
                {
                    indexdata[i] = i;

                    int vi1 = vertexIndices.at(i* 3+ 0 );
                    int uvi1 = uvIndices.at(i* 3+ 0 );
                    int ni1 = normalIndices.at(i* 3+ 0 );

                    Vertex v1;
                    v1.m_position = temp_vertices.at(vi1);
                    v1.m_UV = temp_uvs.at(uvi1);
                    v1.m_normal = temp_normals.at(ni1);

                    int vi2 = vertexIndices.at(i* 3+ 1 );
                    int uvi2 = uvIndices.at(i* 3+ 1 );
                    int ni2 = normalIndices.at(i* 3+ 1 );

                    Vertex v2;
                    v2.m_position = temp_vertices.at(vi2);
                    v2.m_UV = temp_uvs.at(uvi2);
                    v2.m_normal = temp_normals.at(ni2);

                    int vi3 = vertexIndices.at(i * 3 + 2 );
                    int uvi3 = uvIndices.at(i * 3 + 2  );
                    int ni3 = normalIndices.at(i * 3 + 2  );

                    Vertex v3;
                    v3.m_position = temp_vertices.at(vi3);
                    v3.m_UV = temp_uvs.at(uvi3);
                    v3.m_normal = temp_normals.at(ni3);


                    // Calculate the vectors from the current vertex to the two other vertices in the triangle
                    QVector3D v2v1 = v2.m_position - v2.m_position;
                    QVector3D v3v1 = v3.m_position - v1.m_position;

                    // The equation presented in the article states that:

                    // Calculate c2c1_T and c2c1_B
                    float c2c1_T = v2.m_UV.x() - v1.m_UV.x();
                    float c2c1_B = v2.m_UV.y() - v1.m_UV.y();

                    // Calculate c3c1_T and c3c1_B
                    float c3c1_T = v3.m_UV.x() - v1.m_UV.x();
                    float c3c1_B = v3.m_UV.x() - v1.m_UV.x();

                    float fDenominator = c2c1_T * c3c1_B - c3c1_T * c2c1_B;

                    // Calculate the reciprocal value once and for all (to achieve speed)
                    float fScale1 = 1.0f / fDenominator;

                    // T and B are calculated just as the equation in the article states
                    QVector3D  T, B;
                    T = QVector3D((c3c1_B * v2v1.x() - c2c1_B * v3v1.x()) * fScale1,
                            (c3c1_B * v2v1.y() - c2c1_B * v3v1.y()) * fScale1,
                            (c3c1_B * v2v1.z() - c2c1_B * v3v1.z()) * fScale1);

                    B = QVector3D((-c3c1_T * v2v1.x() + c2c1_T * v3v1.x()) * fScale1,
                            (-c3c1_T * v2v1.y() + c2c1_T * v3v1.y()) * fScale1,
                            (-c3c1_T * v2v1.z() + c2c1_T * v3v1.z()) * fScale1);

                    T.normalize();
                    B.normalize();

                    v1.m_tangent = T;
                    v2.m_tangent = T;
                    v3.m_tangent = T;

                    v1.m_bitangent = B;
                    v2.m_bitangent = B;
                    v3.m_bitangent = B;
                }

                g->CreateGeometry(widget,vertexdata,nr_vertex,indexdata,nr_vertex);
                Geometry_List.append(g);
                this->Material_Pointer.append(current_material);

                vertexIndices.clear();
                uvIndices.clear();
                normalIndices.clear();
                unfinished_data = false;
                finish_geometry = false;
            }
        }

        if(set_new_material)
        {
            bool material_set = false;
            for(int i = 0; i < Material_List.length(); i++)
            {
                if(Material_List.at(i)->name.compare(new_material_name) == 0)
                {
                    current_material = i;
                    material_set = true;
                    break;

                }
            }
            if(!material_set)
            {
                current_material = 0;
            }
        }

    }

    if(unfinished_data)
    {
        if(unfinished_data)
        {
            GL3DGeometry * g = new GL3DGeometry();
            int nr_vertex = vertexIndices.size()/3;

            Vertex * vertexdata = (Vertex*) malloc(sizeof(Vertex) * nr_vertex * 3);
            GLuint * indexdata =(GLuint*)  malloc(sizeof(int) * nr_vertex * 3);

            for(int i = 0; i < nr_vertex; i++)
            {
                indexdata[i] = i;

                int vi1 = vertexIndices.at(i* 3+ 0 );
                int uvi1 = uvIndices.at(i* 3+ 0 );
                int ni1 = normalIndices.at(i* 3+ 0 );

                Vertex v1;
                v1.m_position = temp_vertices.at(vi1);
                v1.m_UV = temp_uvs.at(uvi1);
                v1.m_normal = temp_normals.at(ni1);

                int vi2 = vertexIndices.at(i* 3+ 1 );
                int uvi2 = uvIndices.at(i* 3+ 1 );
                int ni2 = normalIndices.at(i* 3+ 1 );

                Vertex v2;
                v2.m_position = temp_vertices.at(vi2);
                v2.m_UV = temp_uvs.at(uvi2);
                v2.m_normal = temp_normals.at(ni2);

                int vi3 = vertexIndices.at(i * 3 + 2 );
                int uvi3 = uvIndices.at(i * 3 + 2  );
                int ni3 = normalIndices.at(i * 3 + 2  );

                Vertex v3;
                v3.m_position = temp_vertices.at(vi3);
                v3.m_UV = temp_uvs.at(uvi3);
                v3.m_normal = temp_normals.at(ni3);


                // Calculate the vectors from the current vertex to the two other vertices in the triangle
                QVector3D v2v1 = v2.m_position - v2.m_position;
                QVector3D v3v1 = v3.m_position - v1.m_position;

                // The equation presented in the article states that:

                // Calculate c2c1_T and c2c1_B
                float c2c1_T = v2.m_UV.x() - v1.m_UV.x();
                float c2c1_B = v2.m_UV.y() - v1.m_UV.y();

                // Calculate c3c1_T and c3c1_B
                float c3c1_T = v3.m_UV.x() - v1.m_UV.x();
                float c3c1_B = v3.m_UV.x() - v1.m_UV.x();

                float fDenominator = c2c1_T * c3c1_B - c3c1_T * c2c1_B;

                // Calculate the reciprocal value once and for all (to achieve speed)
                float fScale1 = 1.0f / fDenominator;

                // T and B are calculated just as the equation in the article states
                QVector3D  T, B;
                T = QVector3D((c3c1_B * v2v1.x() - c2c1_B * v3v1.x()) * fScale1,
                        (c3c1_B * v2v1.y() - c2c1_B * v3v1.y()) * fScale1,
                        (c3c1_B * v2v1.z() - c2c1_B * v3v1.z()) * fScale1);

                B = QVector3D((-c3c1_T * v2v1.x() + c2c1_T * v3v1.x()) * fScale1,
                        (-c3c1_T * v2v1.y() + c2c1_T * v3v1.y()) * fScale1,
                        (-c3c1_T * v2v1.z() + c2c1_T * v3v1.z()) * fScale1);

                T.normalize();
                B.normalize();

                v1.m_tangent = T;
                v2.m_tangent = T;
                v3.m_tangent = T;

                v1.m_bitangent = B;
                v2.m_bitangent = B;
                v3.m_bitangent = B;
            }

            g->CreateGeometry(widget,vertexdata,nr_vertex,indexdata,nr_vertex);
            Geometry_List.append(g);
            this->Material_Pointer.append(current_material);

            vertexIndices.clear();
            uvIndices.clear();
            normalIndices.clear();
            unfinished_data = false;
            finish_geometry = false;
        }

    }

    temp_vertices.clear();
    temp_uvs.clear();
    temp_normals.clear();

    this->Is_Loaded = true;
}

void GL3DModel::Destroy(GL3DWidget * widget)
{

    this->Geometry_List.clear();
    this->Material_List.clear();
    this->Material_Pointer.clear();

    Is_Loaded = false;

}
