Rainfall is modelled from rainfall intensity, which can be spatially and temporally dependent.

$$\frac{\delta V}{\delta t} = P\ s^{2}_{xy} + ...$$

With
$V$ the water volume for a cell $(m^3)$  
$P$ is the rainfall intensity $(m\ h^{-1})$  
$s^{2}_{xy}$ the length and width of the grid cells $(m)$

The cells width and length are corrected for the slope of the local topography. Because of this, on steeper slopes the surface are of the cell becomes larger, and rainfall is spread over a larger area.

$$\delta H = R \cdot \frac{Surface\ length}{cellsize}$$

Snowmelt can be modelled in an identical way, with the snowmelt intensity replacing the rainfall intensity. LISEM does not automatically model snow fall, based on surface temperature. If the user requires to keep track of the depth of the accumulated snow layer, this can be manually done and snowmelt can afterwards be used as input into LISEM.


