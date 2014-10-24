//
//  Utilities.cpp
//  
//
//  Created by Ben Yee on 9/21/13.
//
//

#include "Utilities.h"

vector<double> Utilities::calc_mu_n(int N,double tol){
    vector<double> out;
    double x,x1;
    for(int i = 1; i<=N;i++){
        x = cos(PI*(i-0.25)/(N+0.25));
        do{
            x1 = x;
            x -= lege_eval(N,x)/lege_eval_diff(N,x);
        }while(abs(x1-x)>tol);
        out.push_back(x);
    }
    return out;
}
vector<double> Utilities::calc_w_n(vector<double> mu_n){
    vector<double> out;
    double diff;
    for(int i = 0; i< (int)mu_n.size(); i++){
        diff = lege_eval_diff(mu_n.size(),mu_n[i]);
        out.push_back(2/((1-pow(mu_n[i],2))*pow(diff,2)));
    }
    return out;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
vector<double> Utilities::lege_coef(int N){
    vector< vector<double> > temp;
    for(int i = 0; i<=N; i++){
        vector<double> temp2;
        temp.push_back(temp2);
    }
    //P_0:
    temp[0].push_back(1);
    if(N==0){return temp[0];}
    
    //P_1:
    temp[1].push_back(0);
    temp[1].push_back(1);
    if(N==1){return temp[1];}
    
    //Fill in zeros:
    temp[0].push_back(0);
    for(int i = 2; i<=N; i++){
        temp[0].push_back(0);
        temp[1].push_back(0);
    }
    
    //All the other ones:
    for(int i = 2; i <=N; i++){
        temp[i].push_back(-(i-1)*temp[i-2][0]/i);
        for(int j =1; j<=i;j++){
            temp[i].push_back(((2*i-1)*temp[i-1][j-1]-(i-1)*temp[i-2][j])/i);
        }
        for(int j=i+1;j<=N;j++){
            temp[i].push_back(0);
        }
    }
    return temp[N];
}

double Utilities::lege_eval(vector<double> coeff,double x){
    double out = 0;
    for(int i = 0; i<(int)coeff.size();i++){
        out += coeff[i]*pow(x,i);
    }
    return out;
}


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void Utilities::print_uivector(vector<unsigned int> input_vector,char space){
    cout<<"[";
    for(unsigned int i = 0; i<input_vector.size(); i++){
        cout<<input_vector[i]<<space;
    }
    cout<<"]"<<endl;
}
void Utilities::print_ivector(vector<int> input_vector,char space){
    cout<<"[";
    for(unsigned int i = 0; i<input_vector.size(); i++){
        cout<<input_vector[i]<<space;
    }
    cout<<"]"<<endl;
}
void Utilities::print_dvector(vector<double> input_vector,char space){
    cout<<"[";
    for(unsigned int i = 0; i<input_vector.size(); i++){
        cout<<input_vector[i]<<space;
    }
    cout<<"]"<<endl;
}
void Utilities::print_dmatrix(vector<vector<double> > A, char space){
    cout<<"======Begin matrix output ========"<<endl;
    for(unsigned int i = 0; i<A.size();i++){
        print_dvector(A[i],space);
    }
    cout<<"======End   matrix output ========"<<endl;
}
double Utilities::inf_norm(vector<double> &v1, vector<double> &v2){
    double inf_norm = 0;
    double temp;
    for(unsigned int i = 0; i<min(v1.size(),v2.size()); i++){
        temp = abs(v2[i]-v1[i]);
        if(temp>inf_norm){inf_norm = temp;}
    }
    return inf_norm;
}
double Utilities::p_norm(vector<double> v1, unsigned int p){
    double p_norm = 0;
    unsigned int N = v1.size();
    for(unsigned int i = 0; i<N; i++){
        p_norm += pow(v1[i],p);
    }
    p_norm = p_norm/N;
    p_norm = pow(p_norm,1./p);
    return p_norm;
}
double Utilities::p_norm(vector<double> v1, vector<double> v2, unsigned int p){
    double p_norm = 0;
    unsigned int N = min(v1.size(),v2.size());
    for(unsigned int i = 0; i<N; i++){
        p_norm += pow(v2[i]-v1[i],p);
    }
    p_norm = p_norm/N;
    p_norm = pow(p_norm,1./p);
    return p_norm;
}
double Utilities::p_norm_of_rel_error(vector<double> v_old, vector<double> v_new, unsigned int p, double eps){
    unsigned int N = min(v_old.size(),v_new.size());
    vector<double> rel_error(N,0);
    for(unsigned int i = 0; i < N ; i++ ){
        rel_error[i] = fabs((v_old[i] - v_new[i])/(v_new[i] + eps));
    }
    return Utilities::p_norm(rel_error,p);
}

double Utilities::phi_error(vector<vector<double> > &ref_soln, vector<vector<double> > &soln, int norm){
    unsigned int ref_soln_size = ref_soln[0].size();
    unsigned int soln_size = soln[0].size();
    if(ref_soln_size < soln_size){
        cout<<"!"<<endl;
        cout<<"!"<<endl;
        cout<<"!"<<endl;
        cout<<"ERROR: YOUR REFERENCE SOLUTION MUST BE MORE DETAILED THAN YOUR REGULAR SOLUTION."<<endl;
        cout<<"!"<<endl;
        cout<<"!"<<endl;
        cout<<"!"<<endl;
        return -1;
    }else if(ref_soln[0][0] != soln[0][0] || ref_soln[0][ref_soln_size-1] != soln[0][soln_size-1] ){
        cout<<"!"<<endl;
        cout<<"!"<<endl;
        cout<<"!"<<endl;
        cout<<"ERROR: YOUR REFERENCE SOLUTION AND SOLUTION MUST HAVE MATCHING ENDPOINTS."<<endl;
        cout<<"!"<<endl;
        cout<<"!"<<endl;
        cout<<"!"<<endl;
        return -1;
    }
    vector<double> ref_phi(soln_size,0);
    ref_phi[0] = ref_soln[1][0];
    unsigned int soln_counter = 1;
    for(unsigned int ref_counter = 1; ref_counter<ref_soln_size;ref_counter++){
        if(ref_soln[0][ref_counter] >= soln[0][soln_counter]){
            //Extrapolate...
            ref_phi[soln_counter] = ref_soln[1][ref_counter-1]+
                (ref_soln[1][ref_counter]-ref_soln[1][ref_counter-1])
                *(soln[0][soln_counter]-ref_soln[0][ref_counter-1])/(ref_soln[0][ref_counter]-ref_soln[0][ref_counter-1]);
            soln_counter++;
        }
    }
    
    //Compute p-norm between vectors and return result
    if(norm <= 0){
        return Utilities::inf_norm(ref_phi,soln[1]);
    }else{
        return Utilities::p_norm(ref_phi,soln[1],norm);
    }
}

vector<double> Utilities::vector_add(vector<double> v1, vector<double> v2){
    vector<double> sum;
    for(unsigned int i = 0; i<min(v1.size(),v2.size()); i++){
        sum.push_back(v1[i]+v2[i]);
    }
    return sum;
}
vector<double> Utilities::vector_subtract(vector<double> v1, vector<double> v2){
    vector<double> diff;
    for(unsigned int i = 0; i<min(v1.size(),v2.size()); i++){
        diff.push_back(v1[i]-v2[i]);
    }
    return diff;
}

double Utilities::symmetry_checker(vector<double> &v1){
    double out = 0;
    unsigned int vsize = v1.size();
    for(unsigned int i=0;i<=vsize/2;i++){
        double temp = abs(v1[i]-v1[vsize-i-1]);
        if(temp > out){
            out = temp;
        }
    }
}

bool Utilities::nan_checker(vector<double> &v1){
    for(unsigned int i = 0; i<v1.size();i++){
        if(v1[i]!=v1[i]){
            return true;
        }
    }
    return false;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
vector<double> Utilities::solve_tridiag(vector<vector<double> > A, const vector<double> &b){
    vector<double> x(b);
    unsigned int n = A.size();
    /*
    print_dvector(A[0]);
    print_dvector(A[1]);
    print_dvector(A[2]);
    print_dvector(A[3]);
    print_dvector(A[4]);
    print_dvector(x);
     */
    A[0][2] /= A[0][1];
    x[0] /= A[0][1];
    for(unsigned int i = 1; i<n-1;i++){
        A[i][2] /= A[i][1] - A[i-1][2]*A[i][0];
        x[i] = (x[i]-x[i-1]*A[i][0])/(A[i][1] - A[i-1][2]*A[i][0]);
    }
    unsigned int j = n-1;
    x[j] = (x[j]-x[j-1]*A[j][0])/(A[j][1] - A[j-1][2]*A[j][0]);
    for(int i=n-2; i>=0; i--){
        x[i] -= A[i][2]*x[i+1];
    }
    return x;
}

int Utilities::blowUpChecker(vector<double> &test, double thres){
    for(unsigned int i = 0; i<test.size(); i++){
        if(abs(test[i]) > thres){
            return i;
        }
    }
    return -1;
}

vector<int> Utilities::blowUpChecker(vector<vector<double> > &test, double thres){
    vector<int> out;
    for(unsigned int i = 0; i<test.size(); i++){
        for(unsigned int j = 0; j<test[0].size();j++){
            if(abs(test[i][j]) > thres){
                out.push_back(i);
                out.push_back(j);
                return out;
            }
        }
    }
    out.push_back(-1);
    out.push_back(-1);
    return out;

}

vector<double> Utilities::combine_Phi(const vector<double> &edges, const vector<double> &centers){
    vector<double> combined_vec(edges.size()+centers.size(),0);
    if(edges.size() != centers.size() + 1){
        cout << "WARNING: PROBLEM WITH VECTOR SIZES WHEN COMBINING PHI" << endl;
        return combined_vec;
    }
    combined_vec[0] = edges[0];
    for(unsigned int i = 0; i < centers.size(); i++){
        combined_vec[2*i+1] = centers[i];
        combined_vec[2*i+2] = edges[i+1];
    }
    return combined_vec;
}

void Utilities::split_Phi(const vector<double> &phi_all, vector<double> &phi_edge, vector<double> &phi_cent){
    if(phi_all.size()==0){
        cout<<"Trying to split an empty vector!"<<endl;
        return;
    }
    phi_edge[0] = phi_all[0];
    for(unsigned int i = 0; i<phi_edge.size();i++){
        phi_cent[i] = phi_all[2*i+1];
        phi_edge[i+1] = phi_all[2*(i+1)];
    }
}

double Utilities::find_kappa(double c, const vector<double> &mu_n, const vector<double> &w_n, double tol){
    unsigned int N = w_n.size();
    
    double x = 1.0;
    double x_old = -1000;
    
    while (abs( (x-x_old) / (x+tol*tol) ) > tol){
        x_old = x;
        x -= kappa_fun(c,mu_n,w_n, x)/kappa_fun_deriv(c,mu_n,w_n, x);
    }
    
    return x;
}


double Utilities::kappa_fun(double c, const vector<double> &mu_n, const vector<double> &w_n, double kappa){
    double sum = 2./c;
    for(unsigned int i = 0; i < mu_n.size(); i++){
        sum -= w_n[i]/(1 - kappa*mu_n[i]);
    }
    return sum;
}

double Utilities::kappa_fun_deriv(double c, const vector<double> &mu_n, const vector<double> &w_n, double kappa){
    double sum = 0;
    for(unsigned int i = 0; i < mu_n.size(); i++){
        double term =1 - kappa*mu_n[i];
        sum -= mu_n[i]*w_n[i]/term/term;
    }
    return sum;
}