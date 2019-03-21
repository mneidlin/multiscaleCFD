#include "udf.h"
#include "math.h"


/* Multiscale model of open loop Carotid MS (RC model at ica and eca outlet)
ECA: 0 ICA:1*/

#define nSV 2/*# state variables*/
#define DENSITY 1056.4
real Qeca,Qica;      /*flows from LPN*/
real QFeca,QFica; /*flows calculated by Fluent*/
real Peca,Pica;    /*Pressures from LPN*/
real T,DT,TP1,TP2,INT,STP;

int i,j,NA,NSTEP;


real X[nSV],X1[nSV],X2[nSV]; /*Pressure arrays*/
real DX[nSV],FF[nSV]; /*Delta P and Flow arrays*/
real R[nSV],C[nSV]; /*Resistance and compliance*/



DEFINE_INIT(initialization,domain)
/*remember to put init file (with TP1=TP2=0) and solnet file (with the same left and right columns) in the sim folder*/
{
	#if !RP_NODE	/*SERIAL or HOST*/
		FILE *init,*parameters,;
	    
		
		if ((init=fopen("init.txt","r"))==NULL)
			Message("\nError reading file init.txt\n");
		else
		{
			init=fopen("init.txt","r");
			fscanf(init,"%lf %lf \n %lf %lf %lf  \n",&TP1, &TP2, &INT, &Peca,&Pica);
			fclose(init);
			NSTEP=0;
			Message("initialization \n%lf %lf %lf %lf %lf \n",TP1,TP2,INT,Peca,Pica);
		}
		
/*FROM HERE OWN CODE FOR OPEN LOOP CIRCULATION*/	

/*DECLARATION OF PARAMETERS*/
		if ((parameters=fopen("parameters.txt","r"))==NULL)
			Message("\nError reading file parameters.txt\n");
		else
		{
			parameters=fopen("parameters.txt","r");
			for (i=0;i<nSV;i++)
			fscanf(parameters,"%lf %lf \n",&R[i], &C[i]);
			fclose(parameters);
			
			Message("\n First column: Resistances, Second: Capacitances \n");
			for (i=0;i<nSV;i++)
			Message("%lf %lf\n",R[i], C[i]);
		}
		

	
		
		X2[0]=10; /*Intital pressures eca and ica*/
		X2[1]=10;
			
			

#endif	


}


DEFINE_ADJUST(BOUNDARY,domain)
{
	#if !RP_NODE	/*SERIAL or HOST*/
		int i,FLAG;	
	    real Patr;
	
	#endif
		
	node_to_host_real_2(QFeca,QFica);

	#if !RP_NODE	/*SERIAL or HOST*/
		T = RP_Get_Real("flow-time");
		Message("T: %lf\n",T);		 
	    Patr=30; /*mmHg*/
		
	/*SOLUTION OF ODES WITH REAL TIME*/	
	
	
		if (fabs(T-TP2)<=0.000001)	/*the LPN stops while the 3D model (Fluent) is iterating in the same time step*/
		{
			DT=INT;
			FLAG=2;
		}	
		else if ((T-TP2)>0.000001)	/*jump to the next time step as Fluent has converged*/
		{
			DT=(T-TP2);
			INT=DT;
			TP1=TP2;
			TP2=T;
			FLAG=1;
			NSTEP=NSTEP+1;
		}
		else	/*time step is decreased as Fluent has not converged yet*/
		{
			DT=T-TP1;
			INT=DT;
			TP2=T;
			FLAG=0;
		}
		
		Message("FLAG:%d\n",FLAG);	
	
		if (FLAG!=2)
		{
			Qeca=(QFeca/DENSITY)*pow(10,6);
			Qica=(QFica/DENSITY)*pow(10,6);
		
			
			if (FLAG==1)	/*go to next time-step*/
				for (i=0; i<nSV; i++)
				{
					X[i]=X2[i];	
				}
			else		/*decrease time-step*/
				for (i=0; i<nSV; i++)
				{
					X[i]=X1[i];
				}
				
	/*RC Equations for ECA*/
	
	FF[0]=(X[0]-Patr)/R[0];
	DX[0]=(Qeca-FF[0])/(C[0]);
	
	/*RC Equations for ICA*/
	
	FF[1]=(X[1]-Patr)/R[1];
	DX[1]=(Qica-FF[1])/(C[1]);
	
	
/*Explicit Euler for ODE*/
	X[0]=X[0]+DT*DX[0];
	X[1]=X[1]+DT*DX[1];

	
	
	if (FLAG==1)
	for (i=0; i<nSV; i++)
	{
	X1[i]=X2[i];
	X2[i]=X[i];
	

	}
	else
	for (i=0; i<nSV; i++)
	{
	X2[i]=X[i];
	

	}			
			
			
		
	Peca=(X[0])*133.32;
	Pica=(X[1])*133.32;
	
			
	}	/*end if(FLAG!=2)*/
		
	
		
		
	#endif	/*!RP_NODE*/
	
	host_to_node_real_2(Peca,Pica);

}	/*end DEFINE_ADJUST*/




DEFINE_PROFILE(outlet_eca,thread,nv) /*Pressure BCs, read Fluent MF apply LPN pressure*/
{
	face_t f=0;
	real QF_par=0.0;
	
	 
	#if !RP_HOST	/*SERIAL or NODE*/
		begin_f_loop(f,thread)
			if(PRINCIPAL_FACE_P(f,thread))
			{	
				QF_par += (F_FLUX(f,thread));
				
				F_PROFILE(f,thread,nv) = Peca;
			  
			}
		end_f_loop(f,thread)
	#endif	/*!RP_HOST*/
	
	QFeca=(PRF_GRSUM1(QF_par));
}


DEFINE_PROFILE(outlet_ica,thread,nv) /*Pressure BCs, read Fluent MF apply LPN pressure*/
{
	face_t f=0;
	real QF_par=0.0;
	
	 
	#if !RP_HOST	/*SERIAL or NODE*/
		begin_f_loop(f,thread)
			if(PRINCIPAL_FACE_P(f,thread))
			{	
				QF_par += (F_FLUX(f,thread));
				
				F_PROFILE(f,thread,nv) = Pica;
			  
			}
		end_f_loop(f,thread)
	#endif	/*!RP_HOST*/
	
	QFica=(PRF_GRSUM1(QF_par));
}

DEFINE_PROFILE(mass_flow_cca, thread, position)  /*Pulsatile CCA inflow*/
 {
     
    face_t f;
    real t = CURRENT_TIME;
  real a0 =       6.058 ; 
  real     a1 =     -0.3143  ;
   real    b1 =       2.365  ;
   real    a2 =      -1.095  ;
   real    b2 =       0.905  ;
   real    a3 =     -0.9989  ;
    real   b3 =     -0.3904  ;
    real   a4 =     0.01883  ;
    real   b4 =     -0.4976  ;
    real   a5 =     0.02382  ;
    real   b5 =     -0.4289  ;
    real   a6 =      0.3591  ;
    real   b6 =     0.04926  ;
    real   a7 =     0.06021  ;
    real   b7 =      0.1801  ;
    real   a8 =    -0.08309  ;
    real   b8 =    -0.01789 ; 
	
	real w= 6.219; 
    
	real conversion = 0.0010564;
    begin_f_loop(f,thread)
    {
     
      F_PROFILE(f, thread, position) = conversion*(a0 + a1*cos(t*w) + b1*sin(t*w) + a2*cos(2.*t*w) + b2*sin(2.*t*w) + a3*cos(3.*t*w) + b3*sin(3.*t*w) +  a4*cos(4.*t*w) + b4*sin(4.*t*w) + a5*cos(5.*t*w) + b5*sin(5.*t*w) +a6*cos(6.*t*w) + b6*sin(6.*t*w) + a7*cos(7.*t*w) + b7*sin(7.*t*w) + a8*cos(8.*t*w) + b8*sin(8.*t*w));
    }
    end_f_loop(f, thread)
 
 }

