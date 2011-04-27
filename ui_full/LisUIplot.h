#ifndef LISUIPLOT_H
#define LISUIPLOT_H

#include "lisemqt.h"

//---------------------------------------------------------------------------
class FunctionData: public QwtSyntheticPointData
{
public:
    FunctionData(double(*y)(double)):
        QwtSyntheticPointData(100),
        d_y(y)
    {
    }

    virtual double y(double x) const
    {
        return d_y(x);
    }

private:
    double(*d_y)(double);
};
//---------------------------------------------------------------------------
class Plot : public QwtPlot
{
public:
    Plot( QWidget *parent = NULL);

protected:
    virtual void resizeEvent( QResizeEvent * );

private:
    void populate();
    void updateGradient();
};
//---------------------------------------------------------------------------







#endif // LISUIPLOT_H
