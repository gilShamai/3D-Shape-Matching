/* Compute the Laplace-Beltrami cotangent weight matrix
 * Copyright 2011 Yonathan Aflalo
 */




#include <iostream>
#include <list>
#include <vector>
#include "mex.h"
#include <math.h>
#include "matrix.h"
#include <string.h>
#include <set>

using namespace std;

#define MIN(X, Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X, Y) ((X) > (Y) ? (X) : (Y))
extern int main();

class edge{
public:
    int e1;
    int e2;
    set<int> tri;
    int pos;
    bool operator ==(const edge& b) const{
        return (e1 == b.e1 && e2 == b.e2) || (e1 == b.e2 && e2 == b.e1);
    }
    
    edge(int e1, int e2, int pos) : e1(e1), e2(e2), pos(pos){
    }
    
    edge(int e1, int e2) : e1(e1), e2(e2){
    }
    
    bool operator < (const edge& b) const{
        return  MAX(e1,e2)<MAX(b.e1,b.e2) || (MAX(e1,e2)==MAX(b.e1,b.e2) && MIN(e1,e2)<MIN(b.e1,b.e2));
    }
    int intersect(const edge& b){
        if(e1 == b.e1)
            return e1;
        if(e1 == b.e2)
            return e1;
        if(e2 == b.e1)
            return e2;
        if(e2 == b.e2)
            return e2;
        return -1;
    }
};

class vector3{
public:
    double x;
    double y;
    double z;
    vector3(double x, double y, double z) : x(x), y(y), z(z){
    }
    
    vector3(double x) : x(x), y(x), z(x){
    }
    vector3 operator- () const{
        return vector3(-x, -y, -z);
    }
    friend vector3 operator+ (const vector3 &a, const vector3 &b);
    friend vector3 operator- (const vector3 &a, const vector3 &b);
    friend vector3 operator* (const vector3 &a, const vector3 &b);
    friend vector3 operator* (const double &a, const vector3 &b);
    friend vector3 operator/ (const vector3 &a, const double &b);
    friend double scalar_prod(const vector3 &a, const vector3 &b);
    
    
    double normalize(){
        double norm = sqrt(scalar_prod(*this, *this));
        *this = *this/norm;
        return norm;
    }
    
};

double scalar_prod(const vector3 &a, const vector3 &b){
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

vector3 operator+ (const vector3 &a, const vector3 &b){
    return vector3(a.x + b.x, a.y + b.y, a.z + b.z);
}

vector3 operator- (const vector3 &a, const vector3 &b){
    return vector3(a.x - b.x, a.y - b.y, a.z - b.z);
}
vector3 operator* (const vector3 &a, const vector3 &b){
    return vector3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
}

vector3 operator* (const double &a, const vector3 &b){
    return vector3(a*b.x, a*b.y, a*b.z);
}

vector3 operator/ (const vector3 &a, const double &b){
    return (1/b)*a;
}


const int numInputArgs  = 2;
const int numOutputArgs = 2;


void fill_cotangent_weight(double *sr, mwIndex *irs, mwIndex *jcs, double *tri, int m, size_t size_tri, double *area, double *met);
int sub2ind(int m, int i, int j);
double calcotan(const int &a, const int &b, const int &c, double *met, int i, double *tri);
double calcarea(double *met, int ind_tri, int vertex,double *tri);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != numInputArgs)
        mexErrMsgTxt("Incorrect number of input arguments");
    if (nlhs != numOutputArgs)
        mexErrMsgTxt("Incorrect number of output arguments");
    double* tri_temp = mxGetPr(prhs[0]);
    const size_t size_tri=mxGetM(prhs[0]);
    
    // new: input of G is a matrix of 4 x size_tri
    
    double *G_1D = mxGetPr(prhs[1]);
    double **G = (double **)malloc(size_tri*sizeof(double *));
    
    for(int i = 0; i < size_tri; i++){
        G[i] = &G_1D[i*4];
    }
    
    // previous
//	double **G = new double*[size_tri];
//	for (int i = 0; i< size_tri; i++) {
//		G[i]=new double[4];
//	}
//    for (int i = 0; i< size_tri; i++) {
//        mxArray *G_temp = mxGetFieldByNumber(prhs[1], 1, i);
//        double *mat_temp = mxGetPr(G_temp);
//        for (int j = 0; j<4; j++) {
//            G[i][j]=mat_temp[j];
//        }
//    }
    
    double *met;
	met=new double[3*size_tri];
    mxArray *plhstri[1];
    mxArray *tri_temp_pr[1];
    
    tri_temp_pr[0] = mxCreateDoubleMatrix(size_tri, 3, mxREAL);
    double* tri_temp2 = mxGetPr(tri_temp_pr[0]);
    for(int i=0;i<3*size_tri;i++)
        tri_temp2[i] = tri_temp[i];
    
    
    
    mexCallMATLAB(1, plhstri, 1, &tri_temp_pr[0], "transpose");
    double* tri = mxGetPr(plhstri[0]);
    int m=0;
    for (int i=0; i <size_tri; i++) {
        int a=(int) floor(tri[sub2ind(3, i, 0)])-1;
        int b=(int) floor(tri[sub2ind(3, i, 1)])-1;
        int c=(int) floor(tri[sub2ind(3, i, 2)])-1;
        m = MAX(m,a);
        m = MAX(m,b);
        m = MAX(m,c);
        double lb=sqrt((G[i][3]));
        double lc=sqrt((G[i][0]));
        double la=sqrt((G[i][3]+G[i][0]-2*G[i][2]));
        met[sub2ind(3, i, 0)]=la;
        met[sub2ind(3, i, 1)]=lb;
        met[sub2ind(3, i, 2)]=lc;
    }
    m++;
   
    double maxi_neigh=15;
    mwSize nzmax =(mwSize)ceil((double)m*maxi_neigh);
    mxArray *sprs_mat;
    sprs_mat = mxCreateSparse(m, m, nzmax, mxREAL);
    double *sr  = mxGetPr(sprs_mat);
    mwIndex *irs = mxGetIr(sprs_mat);
    mwIndex *jcs = mxGetJc(sprs_mat);
    double **tri_norm=new double *[(int) size_tri];
    mxArray *area_vec;
    area_vec = mxCreateDoubleMatrix(m, 1, mxREAL);
    double *area=mxGetPr(area_vec);
    fill_cotangent_weight(sr, irs, jcs,tri, m, size_tri, area, met);
    plhs[0] = sprs_mat;
    plhs[1] = area_vec;
    mxDestroyArray(plhstri[0]);
    mxDestroyArray(tri_temp_pr[0]);
    
    free(G);
//    for (int i = 0; i< size_tri; i++) {
//        delete[] G[i];
//    }
//    delete[] G;
    delete[] met;
}

void fill_cotangent_weight(double *sr, mwIndex *irs, mwIndex *jcs, double *tri, const int m, const size_t size_tri, double *area, double *met){
    set<int> *neighboor=new set<int>[m];
    set<int> *triangles=new set<int>[m];
    for(int i=0;i<(int)size_tri;i++){
        int a=(int) floor(tri[sub2ind(3, i, 0)])-1;
        int b=(int) floor(tri[sub2ind(3, i, 1)])-1;
        int c=(int) floor(tri[sub2ind(3, i, 2)])-1;
        neighboor[a].insert(b);
        neighboor[a].insert(c);
        triangles[a].insert(i);
        neighboor[b].insert(a);
        neighboor[b].insert(c);
        triangles[b].insert(i);
        neighboor[c].insert(b);
        neighboor[c].insert(a);
        triangles[c].insert(i);
    }
    int k=0;
    for(int j=0;j<(int) m;j++){
        jcs[j] = k;
        area[j] = 0;
        double sum_cotan=0;
        int index_reach=-1;
        bool flag=1;
        for (set<int>::iterator it=neighboor[j].begin() ; it != neighboor[j].end(); it++ ){
            int current_point=*it;
            set<int> lst_triangle1=triangles[j];
            set<int> lst_triangle2=triangles[current_point];
            set<int> lst_triangle;
            for (set<int>::iterator it2 = lst_triangle1.begin() ; it2 != lst_triangle1.end(); it2++ ){
                set<int>::iterator temp = lst_triangle2.find(*it2);
                if(temp != lst_triangle2.end()){
                    lst_triangle.insert(*temp);
                }
            }
            double curr_cotan=0;
            int sze=0;;
            for (set<int>::iterator it2=lst_triangle.begin() ; it2 != lst_triangle.end(); it2++ ){
                sze++;
                int index_tri = *it2;
                int third_vertex = 0;
                for(int i=0;i<3;i++){
                    if(floor(tri[sub2ind(3, index_tri, i)])-1 != j && floor(tri[sub2ind(3, index_tri, i)])-1 != current_point){
                        third_vertex = (int) floor(tri[sub2ind(3, index_tri, i)])-1;
                        curr_cotan += calcotan(third_vertex, current_point, j, met, index_tri, tri)/2.0;
                        area[j] += calcarea(met, index_tri , j, tri);
                    }
                }
            }
            
            if(flag && current_point>j){
                irs[k] = j;
                sr[k] = 0;
                index_reach = k;
                flag = 0;
                k++;
            }
            sum_cotan+=curr_cotan;
            sr[k] = -curr_cotan;
            irs[k] = current_point;
            k++;
        }
        if(index_reach == -1){
            irs[k] = j;
            sr[k] = 0;
            index_reach = k;
            flag = 0;
            k++;
        }
        sr[index_reach]=sum_cotan;
    }
    jcs[m] = k;
    delete[] triangles;
    delete[] neighboor;
}


void ind2sub(int m, int ind, int *value){
    value[0]=(ind)/m;
    value[1]=ind%m;
}


int sub2ind(int m, int i, int j){
    return i*m+j;
}

double calcotan(const int &a, const int &b, const int &c, double *met, int i, double *tri){
    double l[3];
    int tr[3]={(int)floor(tri[sub2ind(3, i, 0)])-1, (int)floor(tri[sub2ind(3, i, 1)])-1, (int)floor(tri[sub2ind(3, i, 2)])-1};
    
    for (int j = 0; j<3; ++j) {
        if (a == tr[j]) {
            l[2] = met[sub2ind(3, i, j)];
        };
        if (b == tr[j]) {
            l[1] = met[sub2ind(3, i, j)];
        };
        if (c == tr[j]) {
            l[0] = met[sub2ind(3, i, j)];
        };
    }
    double cosalpha = -(l[2] * l[2] - l[0] * l[0] - l[1] * l[1])/(2*l[1]*l[0]);
    
    if (1<fabs(cosalpha)) {
        printf("error\n");
        return 0.5;
    }
    
    double cotalpha=cosalpha/sqrt(1-cosalpha*cosalpha);
    return cotalpha;
}

double calcarea(double *met, int ind_tri, int vertex, double *tri){
    double l1=met[sub2ind(3, ind_tri, 0)];
    double l2=met[sub2ind(3, ind_tri, 1)];
    double l3=met[sub2ind(3, ind_tri, 2)];
    double p = 1/2.0*(l1 + l2 + l3);
    double area_tot = sqrt(p*(p-l1)*(p-l2)*(p-l3));
    double R = (l1 * l2 * l3)/(4.0*area_tot);
    double area = 0;
	int obtuse = -1;
	if (l2*l2+l3*l3<l1*l1)
	{
		obtuse = 0;
	}
	if (l1*l1+l3*l3<l2*l2)
	{
		obtuse = 1;
	}
	if (l1*l1+l2*l2<l3*l3)
	{
		obtuse = 2;
	}
	if (obtuse == -1)
	{
		for (int i = 0; i<3; ++i) {
			if (tri[sub2ind(3, ind_tri, i)]-1 != vertex) {
				double l1_temp = met[sub2ind(3, ind_tri, i)]/2.0;
				double h = sqrt(fabs(R*R - l1_temp * l1_temp)); 
				area += l1_temp * h / 2.0;
			}
		}
	}
	else
	{
		if (tri[sub2ind(3, ind_tri, obtuse)]-1 == vertex)
		{
			area += area_tot/2;
		}
		else
		{
			area += area_tot/4;
		}
	}
	return area;
}
