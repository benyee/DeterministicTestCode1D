
    /*
    if(alpha_mode >=30 || alpha_mode==11){
        if(EDGE_ACCEL_MODE==1){
            Utilities::print_dvector(edgePhi0);
            if(old_phi_0[0] != 0){
                edgePhi0[0] *= phi_0[0]/old_phi_0[0];
            }
            //Update phi_0_edge:
            for(unsigned int i = 1; i < psize;i++){
                if(old_phi_0[i]+old_phi_0[i-1]!=0){
                    edgePhi0[i] *= (phi_0[i]+phi_0[i-1])/(old_phi_0[i]+old_phi_0[i-1]);
                }
            }
            if(old_phi_0[psize-1] != 0){
                edgePhi0[psize] *= phi_0[psize-1]/old_phi_0[psize-1];
            }
        }else if(EDGE_ACCEL_MODE==2){
            edgePhi0[0] = edgePhi0[0] + phi_0[0]-old_phi_0[0];
            //Update phi_0_edge:
            for(unsigned int i = 1; i < psize;i++){
                if(old_phi_0[i]+old_phi_0[i-1]!=0){
                    edgePhi0[i] = edgePhi0[i]+ (phi_0[i]+phi_0[i-1])/2 - (old_phi_0[i]+old_phi_0[i-1])/2;
                }
            }
            edgePhi0[psize] = edgePhi0[psize] + phi_0[psize-1]-old_phi_0[psize-1];
        }else if(EDGE_ACCEL_MODE==3){
            if(old_phi_0[0] != 0){
                edgePhi0[0] *= phi_0[0]/old_phi_0[0];
            }
            //Update phi_0_edge:
            for(unsigned int i = 1; i < psize;i++){
                if(old_phi_0[i]+old_phi_0[i-1]!=0){
                    edgePhi0[i] *= (phi_0[i]/old_phi_0[i]+ phi_0[i-1]/old_phi_0[i-1])/2;
                }
            }
            if(old_phi_0[psize-1] != 0){
                edgePhi0[psize] *= phi_0[psize-1]/old_phi_0[psize-1];
            }
        }else if(EDGE_ACCEL_MODE==4){//HERE I AM
            //This only works if the fine mesh and coarse mesh are the same.
            //This only works if the regions are all uniform
            double lnSth = log(1 + 2.0/sigma_t[0]/h[0]);
            double Sth = sigma_t[0]*h[0];
            double Ssh = sigma_s0[0]*h[0];
            double Aed = (1 - Sth/2*lnSth)/(1 - Ssh/2*lnSth);
            double Bed =3*(Sth*(Sth/2*lnSth-1)+1)/8/(1 - Ssh/2*lnSth);
            double Ced =h[0]/2*lnSth/(1 - Ssh/2*lnSth);
            double Ded =(Sth*(Sth*lnSth-2)+2)/4;
            double Eed =(3*Sth*(Sth*Sth*log(Sth/(2+Sth))+ 2*Sth-2)+8 )/8;
            double Fed =(2 - Sth*lnSth)/4;
            
            double const1 =(Aed-Bed*Ded/Eed)/2;
            double const2 = Bed/Eed;
            
            double denom = 1.0/(Bed*Fed*Ssh/Eed+1);
            
            double Q_mod =(Ced-Bed*Fed*h[0]/Eed)*Q[0];
            
            //First calculate correction factors:
            vector<double> D_cm(edgePhi0);
            if(bc[0] == 1){
                D_cm[0] = (edgePhi0[0] - ( 2*const1*old_phi_0[0] + Q_mod)/denom )/2/old_phi_0[0];
            }else{
                D_cm[0] = (4*in_curr-2*phi_1_e_CM[0])/edgePhi0[0];
            }
            for(unsigned int i = 1; i <psize;i++){
                D_cm[i] = (edgePhi0[i] - ( const1*(old_phi_0[i] + old_phi_0[i-1]) + const2*old_edgePhi1[i] + Q_mod)\
                           /denom )/(old_phi_0[i]+old_phi_0[i-1]);
            }
            if(bc[1] == 1){
                D_cm[psize] = (edgePhi0[psize] - ( 2*const1*old_phi_0[psize-1] + Q_mod)/denom )/2/old_phi_0[psize-1];
            }else{
                D_cm[psize] = (4*out_curr+2*phi_1_e_CM[psize])/edgePhi0[psize];
            }

            Utilities::print_dvector(D_cm);
//            cout<<"Params"<<endl; //HERE I AM
//            cout<<edgePhi0[2]<<endl;
//            cout<<const1*(old_phi_0[2] + old_phi_0[1])/denom<<endl;
//            cout<<const2*old_edgePhi1[2]/denom<<endl;
            
            //Correct:
            if(bc[0] ==1){
                edgePhi0[0] = ( 2*const1*phi_0[0] + Q_mod )/denom + D_cm[0]*phi_0[0]*2;
            }else{
                edgePhi0[0] = (4*in_curr-2*phi_1_e_CM[0])/D_cm[0];
            }
            for(unsigned int i = 1; i <psize;i++){
                edgePhi0[i] = (const1*( phi_0[i] + phi_0[i-1] ) + const2*edgePhi1[i] + Q_mod)\
                    /denom + D_cm[i]*( phi_0[i] + phi_0[i-1] );
            }
            if(bc[1] == 1){
                edgePhi0[psize] = ( 2*const1*phi_0[psize-1] + Q_mod )/denom + D_cm[psize]*phi_0[psize-1]*2;
            }else{
                edgePhi0[psize] = ( 4*out_curr+2*phi_1_e_CM[psize] )/D_cm[psize];
            }
        }
    }*/


    
    if(alpha_mode >=30 || alpha_mode==11){
        if(EDGE_ACCEL_MODE==1){
            unsigned int psize = phi_0.size();
            if(old_phi_0[0] != 0){
                edgePhi0[0] *= phi_0[0]/old_phi_0[0];
                edgePhi1[0] *= phi_0[0]/old_phi_0[0];
            }
            //Update phi_0_edge:
            for(unsigned int i = 1; i < psize;i++){
                if(old_phi_0[i]+old_phi_0[i-1]!=0){
                    edgePhi0[i] *= (phi_0[i]+phi_0[i-1])/(old_phi_0[i]+old_phi_0[i-1]);
                    //This is a cheap way to get edgePhi1.  In reality, you can calculate it from the CMFD equations.
                    edgePhi1[i] *= (phi_0[i]+phi_0[i-1])/(old_phi_0[i]+old_phi_0[i-1]);
                }
            }
            if(old_phi_0[psize-1] != 0){
                edgePhi0[psize] *= phi_0[psize-1]/old_phi_0[psize-1];
                edgePhi1[psize] *= phi_0[psize-1]/old_phi_0[psize-1];
            }
        }else if(EDGE_ACCEL_MODE==2){
            unsigned int psize = phi_0.size();
            edgePhi0[0] = edgePhi0[0] + phi_0[0]-old_phi_0[0];
            if(old_phi_0[0] != 0){
                edgePhi1[0] *= phi_0[0]/old_phi_0[0];
            }
            //Update phi_0_edge:
            for(unsigned int i = 1; i < psize;i++){
                if(old_phi_0[i]+old_phi_0[i-1]!=0){
                    edgePhi0[i] = edgePhi0[i]+ (phi_0[i]+phi_0[i-1])/2 - (old_phi_0[i]+old_phi_0[i-1])/2;
                    //This is a cheap way to get edgePhi1.  In reality, you can calculate it from the CMFD equations.
                    edgePhi1[i] *= (phi_0[i]+phi_0[i-1])/(old_phi_0[i]+old_phi_0[i-1]);
                }
            }
            edgePhi0[psize] = edgePhi0[psize] + phi_0[psize-1]-old_phi_0[psize-1];
            if(old_phi_0[psize-1] != 0){
                edgePhi1[psize] *= phi_0[psize-1]/old_phi_0[psize-1];
            }
        }else if(EDGE_ACCEL_MODE==3){
            unsigned int psize = phi_0.size();
            if(old_phi_0[0] != 0){
                edgePhi0[0] *= phi_0[0]/old_phi_0[0];
                edgePhi1[0] *= phi_0[0]/old_phi_0[0];
            }
            //Update phi_0_edge:
            for(unsigned int i = 1; i < psize;i++){
                if(old_phi_0[i]+old_phi_0[i-1]!=0){
                    edgePhi0[i] *= (phi_0[i]/old_phi_0[i]+ phi_0[i-1]/old_phi_0[i-1])/2;
                    //This is a cheap way to get edgePhi1.  In reality, you can calculate it from the CMFD equations.
                    edgePhi1[i] *= (phi_0[i]+phi_0[i-1])/(old_phi_0[i]+old_phi_0[i-1]);
                }
            }
            if(old_phi_0[psize-1] != 0){
                edgePhi0[psize] *= phi_0[psize-1]/old_phi_0[psize-1];
                edgePhi1[psize] *= phi_0[psize-1]/old_phi_0[psize-1];
            }
        }
    }