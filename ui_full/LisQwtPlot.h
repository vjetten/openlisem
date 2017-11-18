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

    inline void drawCanvas(QPainter * p)
    {

        QwtPlot::drawCanvas(p);

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
