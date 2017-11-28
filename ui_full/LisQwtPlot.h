#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_grid.h>
#include <qwt_plot_marker.h>
#include <qwt_legend.h>
#include <qwt_plot_canvas.h>
#include <qwt_plot_panner.h>
#include <qwt_plot_magnifier.h>
#include <qwt_color_map.h>
#include <qwt_plot_spectrogram.h>
#include <qwt_plot_layout.h>
#include <qwt_matrix_raster_data.h>
#include <qwt_scale_widget.h>
#include <qwt_plot_renderer.h>
#include <qwt_plot_grid.h>
#include <qwt_plot_rescaler.h>
#include <qwt_scale_engine.h>
#include <qwt_plot_zoomer.h>
#include <qwt_picker.h>
#include <qwt_scale_engine.h>
#include <qwt_text.h>
#include <qwt_text_engine.h>
#include <pcrtypes.h>
#include <QMouseEvent>
#include <QPoint>
#include <QWidget>
#include <QList>
#include <QSize>
#include <QLayout>
#include <QVBoxLayout>
#include <QPlainTextEdit>
#include "qwt_plot_curve.h"
#include "qwt_plot_canvas.h"
#include "qwt_plot.h"
#include "global.h"
#include <CsfMap.h>
#include <QMutex>

#ifndef LISQWTPLOT_H
#define LISQWTPLOT_H



class LisQwtPlot : public QwtPlot
{

    public:
    LisQwtPlot( QWidget * w = NULL ) : QwtPlot(w)
    {


    }

    LisQwtPlot( const QwtText &title, QWidget * w= NULL ) : QwtPlot(title,w)
    {


    }

    cTMap * m_velx1;
    cTMap * m_vely1;
    cTMap * m_velx2;
    cTMap * m_vely2;
    double m_startx;
    double m_starty;
    double m_cellsize;
    bool m_velocitySet = false;
    bool m_velocityEnabled = false;
    bool m_velocitySet2 = false;
    double m_velocityPeriod = 25;

    bool profile_enabled = false;
    int profile_currentx = 0;
    int profile_currenty = 0;
    bool profile_started = false;
    bool profile_ended = false;
    bool profile_mousepressed = false;
    QList<QPointF> profile_pointlist;
    QWidget *profile_window= 0;
    QwtPlotCurve *PGraph;
    QwtPlotCurve *MGraph;
    QPlainTextEdit * textGraph;
    QwtPlot * HPlot;

    QMutex map_mutex;
    cTMap * DEM;
    cTMap * DEMChange;

    inline void setVelocityField(cTMap * _velx1, cTMap * _vely1, double cellsize, cTMap * _velx2 = 0, cTMap * _vely2 = 0, double startx = 0, double starty = 0)
    {
        if(_velx1 == 0 || _vely1 == 0)
        {
            return;
        }else
        {
            m_velx1 = _velx1;
            m_vely1 = _vely1;

            m_velocitySet = true;

            m_startx = startx;
            m_starty = starty;
            m_cellsize = cellsize;

            if(_velx2 == 0 || _vely2 == 0)
            {
            }else
            {
                m_velx2 = _velx2;
                m_vely2 = _vely2;

                m_velocitySet2 = true;

            }
        }

    }

    inline void setVelocityFieldEnabled(bool enabled)
    {
        m_velocityEnabled = enabled;

    }

    inline void setVelocityFieldPeriod(int period)
    {
        m_velocityPeriod = period;

    }

    inline void SetProfileMode(bool set)
    {
         profile_enabled = set;

         profile_started = false;
         profile_ended = false;

        profile_pointlist.clear();
    }

    inline void mousePressEvent(QMouseEvent * e)
    {


        QPoint loc = this->canvas()->pos();
        profile_currentx = e->localPos().x() - loc.x();
        profile_currenty = e->localPos().y() - loc.y();



        profile_mousepressed = true;

        replot();

    }

    inline void mouseReleaseEvent(QMouseEvent * e)
    {

        if(e->button() == Qt::MouseButton::LeftButton)
        {
            profile_started = true;

            QPoint loc = this->canvas()->pos();


            QwtInterval interval_x = this->axisInterval(xBottom);
            QwtInterval interval_y = this->axisInterval(yLeft);

            QSize w_size = this->canvas()->size();

            double w_width = w_size.width();
            double w_height = w_size.height();


            profile_pointlist.append(LocalToModelSpace(QPoint(e->localPos().x()-loc.x(),e->localPos().y() - loc.y()),interval_x,interval_y,w_width,w_height));


        }else if(e->button() == Qt::MouseButton::RightButton)
        {
            profile_started = false;
            profile_ended = false;


            StartPlot();

           profile_pointlist.clear();
        }



        profile_mousepressed = false;

        replot();
    }

    inline void StartPlot()
    {

        bool suc = true;
        if(profile_pointlist.length() > 0)
        {

            QList<int> Lx;    //x location
            QList<int> Ly;    //y location
            QList<double> Ls;    //distance from starting point

            QVector<double> SData;
            QVector<double> EData;
            QVector<double> ENData;

            QList<double> Ps;


            for(int i = 0; i < profile_pointlist.length() - 1; i++)
            {
                  double _dx = op.dx;

                    qDebug() << profile_pointlist.at(i).x() << profile_pointlist.at(i+1).x() << "   " << profile_pointlist.at(i).y() <<profile_pointlist.at(i+1).y();

                  double y1 = profile_pointlist.at(i).y() /_dx;
                  double y2 = profile_pointlist.at(i+1).y()/_dx;
                  double x1 = profile_pointlist.at(i).x()/_dx;
                  double x2 = profile_pointlist.at(i+1).x()/_dx;



                  if(y1 == y2 && x1 == x2)
                  {
                      continue;
                  }

                  const bool steep = (fabs(y2 - y1) > fabs(x2 - x1));
                  if(steep)
                  {
                    std::swap(x1, y1);
                    std::swap(x2, y2);
                  }
                  bool reverse = false;
                  if(x1 > x2)
                  {
                    reverse = true;
                    std::swap(x1, x2);
                    std::swap(y1, y2);
                  }

                  const float dx = x2 - x1;
                  const float dy = fabs(y2 - y1);

                  float error = dx / 2.0f;
                  const int ystep = (y1 < y2) ? 1 : -1;
                  int y = (int)y1;

                  const int maxX = (int)x2;

                  for(int x=(int)x1; x<maxX; x++)
                  {
                    if(steep)
                    {

                        Lx.append(y);
                        Ly.append(x);
                        double s = !reverse? std::sqrt(_dx*double((y-y1)*(y-y1)+(x-x1)*(x-x1))) : std::sqrt(_dx*double((y-y2)*(y-y2)+(x-x2)*(x-x2)));
                        for(int k = 0; k < Ps.length(); k++)
                        {
                            s = s +Ps.at(k);
                        }
                        Ls.append(s);
                    }
                    else
                    {

                        Lx.append(x);
                        Ly.append(y);
                        double s = !reverse? std::sqrt(_dx*double((y-y1)*(y-y1)+(x-x1)*(x-x1))) : std::sqrt(_dx*double((y-y2)*(y-y2)+(x-x2)*(x-x2)));
                        for(int k = 0; k < Ps.length(); k++)
                        {
                            s = s +Ps.at(k);
                        }
                        Ls.append(s);
                    }

                    error -= dy;
                    if(error < 0)
                    {
                        y += ystep;
                        error += dx;
                    }
                  }

                  Ps.append(std::sqrt(_dx*double((y2-y1)*(y2-y1)+(x2-x1)*(x2-x1))));

                  map_mutex.lock();
                  for(int j =  0; j < Ls.length() ; j= j+1)
                  {
                      int j2 = reverse? Ls.length() - 1 - j : j;
                      int r = DEM->nrCols() - Ly.at(j2);
                      int c = Lx.at(j2);

                      if(r > 0 && c <DEM->nrCols() && c > 0 && r < DEM->nrRows())
                      {
                          if(!pcr::isMV(op.baseMapDEM->data[r][c]))
                          {
                              SData << Ls.at(j2);
                              EData << DEM->data[r][c];
                              ENData << DEM->data[r][c] + DEMChange->data[r][c];
                          }
                      }
                  }
                  map_mutex.unlock();

                  Ls.clear();
                  Ly.clear();
                  Lx.clear();

            }

            if(SData.length() > 0)
            {




                bool exists = false;
                if(profile_window!= 0)
                {
                    //if(profile_window->isVisible())
                    {
                        exists = true;
                    }
                }

                /*if(exists)
                {

                    //window already exists,
                    //just update the date in the plot


                    PGraph->setSamples(SData,EData);
                    MGraph->setSamples(SData,ENData);

                    textGraph->clear();

                    for(int i = 0; i < Ls.length(); i++)
                    {
                        textGraph->appendPlainText(QString("%1 %2 %3 ")
                                                   .arg(SData.at(i),15,'f',3,' ')
                                                   .arg(EData.at(i),15,'f',3,' ')
                                                   .arg(ENData.at(i),15,'f',3,' '));
                    }

                    HPlot->replot();
                    profile_window->show();

                }else*/
                {
                //create new window


                    profile_window=new QWidget();
                    profile_window->setWindowTitle("profile");
                    profile_window->setMinimumWidth(400);
                    profile_window->setMinimumHeight(500);
                    profile_window->setAttribute(Qt::WA_DeleteOnClose);


                    QwtText title;
                    title.setText("Profile");
                    title.setFont(QFont("MS Shell Dlg 2",12));
                    HPlot = new QwtPlot(title, profile_window);

                    // panning with the left mouse button
                    (void) new QwtPlotPanner( HPlot->canvas() );

                    // zoom in/out with the wheel
                    (void) new QwtPlotMagnifier( HPlot->canvas() );

                    PGraph = new QwtPlotCurve("Elevation");
                    PGraph->attach(HPlot);
                    PGraph->setPen(QPen("#000000"));
                    PGraph->setStyle(QwtPlotCurve::Lines);
                    PGraph->setRenderHint(QwtPlotItem::RenderAntialiased);

                    PGraph->setSamples(SData,EData);

                    MGraph = new QwtPlotCurve("New Elevation");
                    MGraph->attach(HPlot);
                    MGraph->setPen(QPen("#FF0000"));
                    MGraph->setStyle(QwtPlotCurve::Lines);
                    MGraph->setRenderHint(QwtPlotItem::RenderAntialiased);

                    MGraph->setSamples(SData,ENData);


                    HPlot->setCanvasBackground(QBrush(Qt::white));
                    HPlot->enableAxis(HPlot->yRight,true);
                    HPlot->enableAxis(HPlot->yLeft,true);
                    HPlot->setAxisTitle(HPlot->xBottom, "distance (m)");

                    HPlot->setAxisTitle(HPlot->yLeft, "elevation (m)");


                    textGraph = new QPlainTextEdit(this);
                    textGraph->setWordWrapMode(QTextOption::NoWrap);
                    textGraph->setMaximumHeight(80);
                    textGraph->clear();

                    for(int i = 0; i < SData.length(); i++)
                    {
                        textGraph->appendPlainText(QString("%1 %2 %3 ")
                                                   .arg(SData.at(i),15,'f',3,' ')
                                                   .arg(EData.at(i),15,'f',3,' ')
                                                   .arg(ENData.at(i),15,'f',3,' '));
                    }
                    QVBoxLayout *layout= new QVBoxLayout(profile_window);
                    layout->addWidget(HPlot);
                    layout->addWidget(textGraph);

                    profile_window->show();
                    HPlot->replot();
                }

        }else
        {
            suc = false;
        }


        }else
        {
            suc = false;
        }

        if(!suc)
        {

            QMessageBox::warning(this, "openLISEM",
                                 QString("not enough points in map extent for profile")
                                 );
        }
    }

    inline QPointF LocalToModelSpace(QPoint p,QwtInterval i_x,QwtInterval i_y, int width, int height)
    {
        QPointF ret;

        ret.setX(i_x.minValue() + (double(p.x())/double(width)) * (i_x.maxValue() - i_x.minValue()));
        ret.setY(i_y.minValue() + ((double(height) -double(p.y()))/double(height)) * (i_y.maxValue() - i_y.minValue()));

        return ret;
    }

    inline QPointF ModelToLocalSpace(QPointF p,QwtInterval i_x,QwtInterval i_y, int width, int height)
    {
        QPointF ret;

        ret.setX(double(width) * (double(p.x() - i_x.minValue()))/(i_x.maxValue() - i_x.minValue()));
        ret.setY(double(height) - double(height) * (double(p.y() - i_y.minValue()))/(i_y.maxValue() - i_y.minValue()));

        return ret;

    }


    inline void mouseMoveEvent(QMouseEvent * e)
    {

        QPoint loc = this->canvas()->pos();

        profile_currentx = e->localPos().x()-loc.x();
        profile_currenty = e->localPos().y() - loc.y();


        replot();

    }

    inline void mouseDoubleClickEvent(QMouseEvent * e)
    {


        replot();
    }

    inline void drawItems(QPainter *p, const QRectF & r,
                          const QwtScaleMapTable & s)
    {
        QwtPlot::drawItems(p,r,s);


    };

    inline void drawCanvas(QPainter * p)
    {

        //draw original canvas first (the map display as original
        QwtPlot::drawCanvas(p);



        //draw Lines for the profile tool
        //draw from all existing lines in the point list
        //finally, if currently picking another point, draw towards the mouse cursor
        if(profile_enabled && profile_started == true && profile_ended == false)
        {

            QwtInterval interval_x = this->axisInterval(xBottom);
            QwtInterval interval_y = this->axisInterval(yLeft);

            QSize w_size = this->canvas()->size();

            double w_width = w_size.width();
            double w_height = w_size.height();

            for(int i = 0; i < profile_pointlist.length() - 1; i++)
            {
                QPointF p1 = ModelToLocalSpace(profile_pointlist.at(i),interval_x,interval_y,w_width,w_height);
                QPointF p2 = ModelToLocalSpace(profile_pointlist.at(i+1),interval_x,interval_y,w_width,w_height);

                p->setPen(QColor(255,0,0,255));
                p->drawLine(p1,p2);

            }

            if(profile_pointlist.length()>0 && profile_mousepressed)
            {
                QPointF p1 = ModelToLocalSpace(profile_pointlist.at(profile_pointlist.length() - 1),interval_x,interval_y,w_width,w_height);
                QPointF p2 = QPointF(profile_currentx,profile_currenty);

                p->setPen(QColor(255,0,0,255));
                p->drawLine(p1,p2);

            }

        }



        //Draw velocity field
        //
        //basicly iterate through the canvas, for every nth pixel (determined by velocity_period), draw a line in the direction of the velocity map
        if(m_velocityEnabled && m_velocitySet && m_velx1 != 0 && m_vely1 != 0)
        {
            double max_velsqr = 0;

            for(int i = 0; i < m_velx1->nrCols() ; i++)
            {
                for(int j = 0; j < m_velx1->nrRows() ; j++)
                {

                    if(!pcr::isMV(m_velx1->data[j][i]))
                    {
                        max_velsqr = std::max(max_velsqr, m_velx1->data[j][i] *m_velx1->data[j][i] + m_vely1->data[j][i] *m_vely1->data[j][i]);
                    }
                }
            }
            double max_vel = std::sqrt(max_velsqr);

            QwtInterval interval_x = this->axisInterval(xBottom);
            QwtInterval interval_y = this->axisInterval(yLeft);

            QSize w_size = this->canvas()->size();

            double w_width = w_size.width();
            double w_height = w_size.height();

            double p_endx = interval_x.maxValue();
            double p_endy = interval_y.maxValue();

            double p_startx = interval_x.minValue();
            double p_starty = interval_y.minValue();

            int n_x = std::floor(w_width/m_velocityPeriod);
            int n_y = std::floor(w_height/m_velocityPeriod);

            double drawstartx = 0;
            double drawstarty = 0;

            for(int i = 0; i < n_x; i++)
            {
                for(int j = 0; j < n_y; j++)
                {
                    double l_x = i * m_velocityPeriod;
                    double l_y = j * m_velocityPeriod;

                    double p_x = p_startx + (p_endx - p_startx)*(l_x/w_width);
                    double p_y = m_velx1->nrRows() * m_cellsize - (p_endy - (p_endy - p_starty)*(l_y/w_height));

                    QPointF vel = getVelocityAt(m_velx1,m_vely1,p_x,p_y);
                    double vel_l = std::sqrt(double(vel.x()*vel.x() + vel.y()*vel.y()));

                    if(vel_l > 0.001)
                    {
                        double d_x =  2.0*(double(vel.x())/vel_l)* (vel_l/max_vel)* m_velocityPeriod*2.0/3.0;
                        double d_y = 2.0*(double(vel.y())/vel_l)* (vel_l/max_vel) * m_velocityPeriod*2.0/3.0;

                        if(std::fabs(d_x) > 1.0 || std::fabs(d_y) > 1.0)
                        {
                            p->setPen(QColor(0,0,0,255));
                            p->drawLine(QPoint(l_x,l_y),QPoint(l_x + d_x, l_y + d_y));
                        }
                    }


                }
            }

        }
    }

    inline QPointF getVelocityAt(cTMap * _velx, cTMap * _vely,double x, double y)
    {

        double cs = _velx->cellSize();
        int cn = std::floor(x/cs);
        int rn = std::floor(y/cs);

        double v1 = 0;
        double v2 = 0;
        double v3 = 0;
        double v4 = 0;

        double u1 = 0;
        double u2 = 0;
        double u3 = 0;
        double u4 = 0;

        double w1 = 0;
        double w2 = 0;
        double w3 = 0;
        double w4 = 0;

        if(!OUTOFMAP(_velx,rn,cn))
        {
            if(!MV(_velx,rn,cn))
            {
                w1 = 1;
                v1 =_vely->data[rn][cn];
                u1 =_velx->data[rn][cn];
            }
        }
        if(!OUTOFMAP(_velx,rn+1,cn))
        {
            if(!MV(_velx,rn+1,cn))
            {
                w2 = 1;
                v2 =_vely->data[rn+1][cn];
                u2 =_velx->data[rn+1][cn];
            }
        }
        if(!OUTOFMAP(_velx,rn,cn+1))
        {
            if(!MV(_velx,rn,cn+1))
            {
                w3 = 1;
                v3 =_vely->data[rn][cn+1];
                u3 =_velx->data[rn][cn+1];
            }
        }
        if(!OUTOFMAP(_velx,rn+1,cn+1))
        {
            if(!MV(_velx,rn+1,cn+1))
            {
                w4 = 1;
                v4 =_vely->data[rn+1][cn+1];
                u4 =_velx->data[rn+1][cn+1];
            }
        }

        QPointF vel = QPoint(0,0);
        vel.setX((w1+w2+w3+w4) > 0? ((w1*u1 + w2*u2 + w3*u3 + w4 * u4)/(w1+w2+w3+w4)): 0.0);
        vel.setY((w1+w2+w3+w4) > 0? ((w1*v1 + w2*v2 + w3*v3 + w4 * v4)/(w1+w2+w3+w4)): 0.0);

        return vel;
    }

    inline void drawArrow(QPainter * p, double x_or, double y_or, double x_dir, double y_dir, double length)
    {




    }



};


#endif // LISQWTPLOT_H
