#include "sode.h"

sode::sode(sode_type type, int dim){
    this->type = type;
    this->dim = dim;
	this->verb = false;

	this->bisecting = false;
	this->h = 1e-3;
    this->h_min = 1e-6;
    this->h_max = 1e-1;
	this->tol_end = 1e-15;
    this->tol_sol = 1e-10; 
	this->tol_mon = 1e-12;
    this->damp = 0.9;

    this->x0 = new double[dim];
	this->x1 = new double[dim];
	this->x2 = new double[dim];

    switch (type){
    case SODE_RK56_FB: // Fehlberg
        alloc_56FB();
        break;
    case SODE_RK56_CK: // Cash-Karp
        alloc_56CK();
        break;
    case SODE_RK78_DP: // Dormand-Prince
        alloc_RK78();
        break;
    default:
        break;
    }
}

sode::~sode(){
	delete(x0);
	delete(x1);
	delete(x2);
	delete(b1);
	delete(b2);
	delete(c);

	for (int i=0; i<s; i++){
		 delete(a[i]);
		 delete(k[i]);
	}
	delete(a);
	delete(k);
}

void
sode::configure(double h_init, double h_min, double h_max){
	this->h_init = h_init;
    this->h_min = h_min;
    this->h_max = h_max;
	this->h = h_min;
}
void
sode::configure(double tol_sol, double tol_mon, double tol_end, double damp){
	this->tol_sol = tol_sol;
	this->tol_mon = tol_mon;
	this->tol_end = tol_end;
	this->damp = damp;
}

void
sode::set_system(int (*system)(double* f, double* x, double t, void* aux)){
    this->system = system;
}

void
sode::set_monitor(bool(*monitor)(double* x, double t, void* aux)){
    this->monitor = monitor;
}

void
sode::reset(){
	this->h = h_min;
	this->bisecting = false;
}

int
sode::evolve(double* x, double* t, double t_end, int chg_mon, void* aux){
	if (verb) printf("sode::evolve\n");
    int status = calc_ks (x, *t, aux);
    if (status != SODE_OK){
        h = 0.5*h;
		if (verb) printf("bad function at t = %f\n", *t);
        return status; // return SODE_BAD_FUNC
    }
	if (verb) printf("sode::evolve: 1\n");
	// calculate high and low order steps
    vector_copy (x1, x, dim);
    vector_copy (x2, x, dim);
    for (int i=0; i<s; i++){
        vector_combine (x1, 1.0, h*b1[i], x1, k[i], dim);
	    vector_combine (x2, 1.0, h*b2[i], x2, k[i], dim);
    }
	if (verb) printf("sode::evolve: 2\n");
    double err = vector_dist_chbyv (x1, x2, dim);
    if (err > tol_sol){
		if (h < h_min){
			if (verb) printf("stepsize falled below hmin: %.3e\n", h_min);
			return SODE_FAILED;
		}else {
        	h = 0.5*h;
			if (verb) printf("err = %.3e, reduced step-size to %.3e\n", err, h);
        	return SODE_CONTINUE_BAD_STEP;
		}
    }
	if (verb) printf("sode::evolve: 3\n");
    // event detection and bisection
    if (chg_mon != 0 && monitor){
		if (verb) printf("sode::evolve: 30\n");
        bool mon0 = monitor(x, *t, aux);
        bool mon1 = monitor(x1, *t + h, aux);
        if (mon1 - mon0 == chg_mon){
			if (verb) printf("sode::evolve: 31\n");
			bisecting = true;
            // check tolerance
			double eps = vector_dist_chbyv (x, x1, dim);
            if (eps > tol_mon){
				if (verb) printf("sode::evolve: 32\n");
                // bisect stepsize and return
                h = 0.5*h;
				if (verb) printf("sode::evolve: eps = %f, bisected h: %.3e\n", eps, h);
                return SODE_CONTINUE_BAD_STEP;
            } else {
				if (verb) printf("sode::evolve: 33\n");
				bisecting = false;
				printf("sode::evolve: event found, eps = %.2e\n", eps);
                vector_copy (x, x1, dim);
                *t = *t + h;
                return SODE_SUCCESS_MONITOR;
            }
        }
    }
	if (verb) printf("sode::evolve: 4\n");
    // perform step
    vector_copy (x, x1, dim);
    *t = *t + h;

    if (!bisecting){// regular step-size control
    	double r;
    	if (err > 0.5*tol_sol) r = damp * pow (tol_sol/err, 1.0/q);
		else r = damp * pow (tol_sol/err, 1.0/(q+1.0));
		// bound r, h: prevent feedback oscillations & halting
		if (r < 0.2) r = 0.2;
		if (r > 5.0) r = 5.0;
		h =  r * h;
		if (fabs(h) < h_min) h = h_min*fabs(h)/h;
    	if (fabs(h) > h_max) h = h_max*fabs(h)/h;
		if (verb) printf("regular step control h: %.3f, r: %.3e\n", h, r);
	}
	if (verb) printf("sode::evolve: 5\n");
    if (fabs(*t - t_end) < tol_end) return SODE_SUCCESS_TIME;
	else if (fabs(*t - t_end) < h) h = t_end - *t;

    return SODE_CONTINUE_GOOD_STEP;
}

void
sode::set_verb(){
	this->verb = true;
}


int
sode::calc_ks(double* x, double t, void* aux){
    int status = 0;
    for (int i=0; i<s; i++){
		vector_copy(x0, x, dim);
	    for (int j=0; j<i; j++){
		    vector_combine (x0, 1.0, h*a[i][j], x0, k[j], dim);
	    }
	    status = system (k[i], x0, t+c[i]*h, aux);
        if (status != 0) break;
	}
    if (status != 0) return SODE_BAD_FUNC;
    else return SODE_OK;
}

void
sode::vector_combine(double* x_comb, double a, double b,
				     double* xa, double* xb, int size){
    for (int i=0; i<size; i++) x_comb[i] = a*xa[i] + b*xb[i];
}

void
sode::vector_copy (double* x_copy, double* x, int size){
    for (int i=0; i<size; i++) x_copy[i] = x[i];
}

double
sode::vector_dist_chbyv (double* xa, double* xb, int size){
    double dist = 0;
    for (int i=0; i<size; i++){
        if (fabs(xa[i]-xb[i]) > dist) dist = fabs(xa[i]-xb[i]);
    }
    return dist;
}

void
sode::alloc_56FB(){
    q = 5;
	s = 6;
	
	c = new double[s];
	c[0]=0.0;     c[1]=1./4.;  c[2]=3./8.;
	c[3]=12./13.; c[4]=1.;     c[5]=1./2.;
	
	b1 = new double[s];
	b1[0]=16./135.;      b1[1]=0.;      b1[2]=6656./12825.;
	b1[3]=28561./56430.; b1[4]=-9./50.; b1[5]=2./55.;
	
	b2 = new double[s];
	b2[0]=25./216.;      b2[1]=0.;      b2[2]=1408./2565.;
	b2[3]=2197./4104.;   b2[4]=-1./5.;  b2[5]=0.;
	
	a = new double*[s];
	for (int i=0; i<6; i++) a[i] = new double[i];
	
	a[1][0]=1./4.;
	a[2][0]=3./32.;      a[2][1]=9./32.;
	
	a[3][0]=1932./2197.; a[3][1]=-7200./2197.;
	a[3][2]=7296./2197.;
	
	a[4][0]=439./216.;   a[4][1]=-8.;
	a[4][2]=3680./513.;  a[4][3]=-845./4104.;
	
	a[5][0]=-8./27.;     a[5][1]=2.;        a[5][2]=-3544./2565.;
	a[5][3]=1859./4104.; a[5][4]=-11./40.;
        
	k = new double*[s];
	for (int i=0; i<s; i++) k[i] = new double[dim];
}

void
sode::alloc_56CK(){
    q = 5;
	s = 6;
	
	c = new double[s];
	c[0]=0.0;     c[1]=1./5.;  c[2]=3./10.;
	c[3]=3./5.;   c[4]=1.;     c[5]=7./8.;
	
	b1 = new double[s];
	b1[0]=37./378.;      b1[1]=0.;      b1[2]=250./621.;
	b1[3]=125./594.;     b1[4]=0.;      b1[5]=512./1771.;
	
	b2 = new double[s];
	b2[0]=2825./27648.;  b2[1]=0.;          b2[2]=18575./48384.;
	b2[3]=13525./55296.; b2[4]=277./14336.; b2[5]=1./4.;
	
	a = new double*[s];
	for (int i=0; i<6; i++) a[i] = new double[i];
	
	a[1][0]=1./5.;
	a[2][0]=3./40.;          a[2][1]=9./40.;
	
	a[3][0]=3./10.;          a[3][1]=-9./10.;     a[3][2]=6./5.;
	
	a[4][0]=-11./54.;        a[4][1]=5/2.;        a[4][2]=-70./27.;
	a[4][3]=35./27.;
	
	a[5][0]=1631./55296.;   a[5][1]=175./512.;  a[5][2]=575./13824.;
	a[5][3]=44275./110592.; a[5][4]=253./4096.;
        
	k = new double*[s];
	for (int i=0; i<s; i++) k[i] = new double[dim];
}

void
sode::alloc_RK78(){
    q = 7;
	s = 13;
	
	c = new double[s];
	c[0] = 0.;	         	        	c[1] = 1./18.;
	c[2] = 1./12.;		         		c[3] = 1./8.;
	c[4] = 5./16.;		       	 		c[5] = 3./8.;
	c[6] = 59./400.;		        	c[7] = 93./200.;
	c[8] = 5490023248./9719169821.;		c[9] = 13./20.;
	c[10] = 1201146811./1299019798.;	c[11] = 1.;
	c[12] = 1.;
	
	b1 = new double[s];
	b1[0] = 14005451./335480064.;    b1[1] = 0.;
	b1[2] = 0.;         		    b1[3] = 0.;
	b1[4] = 0.;       		    b1[5] = -59238493./1068277825.;
	b1[6] = 181606767./758867731.;   b1[7] = 561292985./797845732.;
	b1[8] = -1041891430./1371343529.;b1[9] = 760417239./1151165299.;
	b1[10] = 118820643./751138087.;  b1[11] = -528747749./2220607170;
	b1[12] = 0.25;
	
	b2 = new double[s];
	b2[0] = 13451932./455176623.;       b2[1] = 0.;
	b2[2] = 0.;       		       b2[3] = 0.;
	b2[4] = 0.;      		       b2[5] = -808719846./976000145.;
	b2[6] = 1757004468./5645159321.;    b2[7] = 656045339./265891186.;
	b2[8] = -3867574721./1518517206.;   b2[9] = 465885868./322736535.;
	b2[10] = 53011238./667516719.;      b2[11] = 2./45.;
	b2[12] = 0.;
	
	a = new double*[s];
	for (int i=0; i<13; i++) a[i] = new double[i];
	a[1][0] = 1./18.;
	a[2][0] = 1./48.;
	a[3][0] = 1./32.;
	a[4][0] = 5./16.;
	a[5][0] = 3./80.;
	a[6][0] = 29443841./614563906.;
	a[7][0] = 16016141./946692911.;
	a[8][0] = 39632708./573591083.;
	a[9][0] = 246121993./1340847787.;
	a[10][0] =	-1028468189./846180014.;
	a[11][0] = 185892177./718116043.;
	a[12][0] = 403863854./491063109.;
	
	a[2][1] = 1./16.;
	a[3][1] = 0.;
	a[4][1] = 0.;
	a[5][1] = 0.;
	a[6][1] = 0.;
	a[7][1] = 0.;
	a[8][1] = 0.;
	a[9][1] = 0.;
	a[10][1] = 0.;
	a[11][1] = 0.;
	a[12][1] = 0.;
	
	a[3][2] = 3./32.;
	a[4][2] = -75./64.;
	a[5][2] = 0.;
	a[6][2] = 0.;
	a[7][2] = 0.;
	a[8][2] = 0.;
	a[9][2] = 0.;
	a[10][2] = 0.;
	a[11][2] = 0.;
	a[12][2] = 0.;
	
	a[4][3] = 75./64.;
	a[5][3] = 3./16.;
	a[6][3] = 77736538./692538347.;
	a[7][3] = 61564180./158732637.;
	a[8][3] = -433636366./683701615.;
	a[9][3] = -37695042795./15268766246.;
	a[10][3] = 8478235783./508512852.;
	a[11][3] = -3185094517./667107341.;
	a[12][3] = -5068492393./434740067.;
	
	a[5][4] = 3./20.;
	a[6][4] = -28693883./1125000000.;
	a[7][4] = 22789713./633445777.;
	a[8][4] = -421739975./2616292301.;
	a[9][4] = -309121744./1061227803.;
	a[10][4] = 1311729495./1432422823.;
	a[11][4] = -477755414./1098053517.;
	a[12][4] = -411421997./543043805.;
	
	a[6][5] = 23124283./1800000000.;
	a[7][5] = 545815736./2771057229.;
	a[8][5] = 100302831./723423059.;
	a[9][5] = -12992083./490766935.;
	a[10][5] = -10304129995./1701304382.;
	a[11][5] = -703635378./230739211.;
	a[12][5] = 652783627./914296604.;
	
	a[7][6] = -180193667./1043307555.;
	a[8][6] = 790204164./839813087.;
	a[9][6] = 6005943493./2108947869.;
	a[10][6] = -48777925059./3047939560.;
	a[11][6] = 5731566787./1027545527.;
	a[12][6] = 11173962825./925320556.;
	
	a[8][7] = 800635310./3783071287.;
	a[9][7] = 393006217./1396673457.;
	a[10][7] = 15336726248./1032824649.;
	a[11][7] = 5232866602./850066563.;
	a[12][7] = -13158990841./6184727034.;
	
	a[9][8] = 123872331./1001029789.;
	a[10][8] = -45442868181./3398467696.;
	a[11][8] = -4093664535./808688257.;
	a[12][8] = 3936647629./1978049680.;
	
	a[10][9] = 3065993473./597172653.;
	a[11][9] = 3962137247./1805957418.;
	a[12][9] = -160528059./685178525.;
	
	a[11][10] = 65686358./487910083.;
	a[12][10] = 248638103./1413531060.;
	
	a[12][11] = 0.;
	
	k = new double*[s];
	for (int i=0; i<s; i++) k[i] = new double[dim];
}
