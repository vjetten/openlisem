/* headcalc.h                                    */
/*  compute theta, k and dif moist cap from head */



double TheNode(
		double head,
		const  HORIZON *hor);

double HcoNode(
		double head,
		const HORIZON *hor,
		double calib);

double DmcNode(
		double head,
		const  HORIZON *hor);

