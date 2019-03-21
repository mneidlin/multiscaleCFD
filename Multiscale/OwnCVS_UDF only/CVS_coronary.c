#include "udf.h"
#include "math.h"


#define nSV 13	/*# state variables*/
#define DENSITY 1056.4



real Qlver,Qlsub,Qlcar,Qrver,Qrsub,Qrcar,Qlcor1,Qlcor2,Qrcor,Qao,Qheart;      /*flows from LPN*/
real QFlver,QFlsub,QFlcar,QFrver,QFrsub,QFrcar,QFlcor1,QFlcor2,QFrcor,QFao,QFheart; /*flows calculated by Fluent*/
real Pub,Plb,Pheart,Pcor;       /*Pressures from LPN*/
real T,DT,TP1,TP2,INT,STP;

int i,j,NA,NSTEP;

/*From here OWN Variables*/
/*
TIMESTEP OF 1e-4 is necessary for Explicit! TS of 2e-3 up to 5e-3 is possible for implicit/explicit mix !!! 
%values of R (mmHg ml/s) and C ml/mmHg and CQi (ml /s mmHg^0.5) and Volumes
%of Compartments V0 (ml)
%1: sart_UB 2: svn_UB 3: sart_lb, 4:svn_LB, 5:RA, 6:RV, 7:pulm_artery 8:
%part 9: pvn 10: LA 11:LV 12:COR 

*/

real V[nSV],V1[nSV],V2[nSV]; /*Volume arrays (Heart)*/
real X[nSV],X1[nSV],X2[nSV]; /*Pressure arrays, X(5,6,10,11) are IRRELEVANT , Ventricles Atria computed with V array!*/ 
real FF[nSV], PP[nSV]; /*Flows through resistances and Pressures at ventricles7atria*/
real DV[nSV],DX[nSV];
real R[nSV],C[nSV]; /*Resistance and compliance*/
real EmaxRA,EdRA,EmaxRV,EmaxLA,EdLA,EmaxLV,LV_Pd_beta,LV_Pd_kappa,LV_Pd_alpha,RV_Pd_beta,RV_Pd_kappa,RV_Pd_alpha; /*Elastances of Ventricles and atria*/


DEFINE_INIT(initialization,domain)
/*remember to put init file (with TP1=TP2=0) and solnet file (with the same left and right columns) in the sim folder*/
{
	#if !RP_NODE	/*SERIAL or HOST*/
		FILE *init,*solnet,*parameters,*heart;
	    
		
		if ((init=fopen("init.txt","r"))==NULL)
			Message("\nError reading file init.txt\n");
		else
		{
			init=fopen("init.txt","r");
			fscanf(init,"%lf %lf \n %lf %lf %lf %lf\n",&TP1, &TP2, &INT, &Pub, &Plb, &Pheart);
			fclose(init);
			NSTEP=0;
			Message("initialization \n%lf %lf %lf %lf %lf %lf\n",TP1,TP2,INT,Pub,Plb,Pheart);
		}
		
/*FROM HERE OWN CODE FOR CLOSED LOOP CIRCULATION*/	

/*DECLARATION OF PARAMETERS*/
		if ((parameters=fopen("parameters.txt","r"))==NULL)
			Message("\nError reading file parameters.txt\n");
		else
		{
			parameters=fopen("parameters.txt","r");
			for (i=1;i<nSV;i++)
			fscanf(parameters,"%lf %lf \n",&R[i], &C[i]);
			fclose(parameters);
			
			Message("\n First WORSCHT61 column: Resistances, Second: Capacitances \n");
			for (i=1;i<nSV;i++)
			Message("%lf %lf\n",R[i], C[i]);
		}
		
/*Heart Parameters*/		
		if ((heart=fopen("heart.txt","r"))==NULL)
			Message("\nError reading file heart.txt\n");
		else
		{
			heart=fopen("heart.txt","r");
			fscanf(heart,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&EmaxRA, &EdRA, &EmaxRV, &EmaxLA, &EdLA, &EmaxLV,&LV_Pd_beta,&LV_Pd_kappa,&LV_Pd_alpha,&RV_Pd_beta,&RV_Pd_kappa,&RV_Pd_alpha);
			fclose(heart);
			
			Message("\n EmaxRA,EdRA,EmaxRV,EmaxLA,EdLA,EmaxLV LV_Pd_beta,LV_Pd_kappa,LV_Pd_alpha,RV_Pd_beta,RV_Pd_kappa,RV_Pd_alpha\n");
			Message("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",EmaxRA, EdRA, EmaxRV, EmaxLA, EdLA, EmaxLV,LV_Pd_beta,LV_Pd_kappa,LV_Pd_alpha,RV_Pd_beta,RV_Pd_kappa,RV_Pd_alpha);
		}
		
/*Needed for restart of simulation during crash*/
		if ((solnet=fopen("solnet_values.txt","r"))==NULL) 
			Message("\n Error reading file solnet_values.txt\n");
		else
		{
			solnet=fopen("solnet_values.txt","r");
			for (i=1;i<nSV;i++)
			fscanf(solnet,"%lf %lf %lf \n",&V[i], &V2[i], &V1[i]);
			fclose(solnet);
			
			Message("\n Initial Volumes \n");
			for (i=1;i<nSV;i++)
			Message("%lf \n",V2[i]);
			
			for (i=1;i<nSV;i++)
			X2[i]=V2[i]/C[i];
		
			X2[0]=100;
			
			
		Message("\n Initial Pressures \n");
			for (i=1;i<nSV;i++)
			Message("%lf \n",X2[i]);
		}
		

#endif	


}


DEFINE_ADJUST(BOUNDARY,domain)
{
	#if !RP_NODE	/*SERIAL or HOST*/
		FILE *solnet,*results,*pv;
		real modulo;
		real ev,ea,Ts1,Ts2,Tsa1,Tsa2,Zc1,Zc2;
		int i,FLAG;	
		real Qub,Qlb,Qcor;
		Ts1=0.26;
		Ts2=0.39;
		Tsa1=0.04;
		Tsa2=0.09;
		
		Zcor=1.4; /*Theoretically it should be 1.4 to allow a 95%, 5% Flow split!*/
		Zub=0.2; /*Additional resistance representing the impedance of the 3D model, Z=c*rho/A c:PWV=sqrt(Vaorta/(rho*Caorta)), was additionally adjusted with Matlab CVS! 0.2 instead of 0.28*/
		Zlb=0.112;
	
	#endif
		
	node_to_host_real_11(QFlver,QFlsub,QFlcar,QFrver,QFrsub,QFrcar,QFlcor1,QFlcor2,QFrcor,QFao,QFheart);
		
	#if !RP_NODE	/*SERIAL or HOST*/
		T = RP_Get_Real("flow-time");
		Message("T: %lf\n",T);		 
	
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
			
			Qlver=(QFlver/DENSITY)*pow(10,6); /*conversion from kg/s to ml/s */
			Qlcar=(QFlcar/DENSITY)*pow(10,6);
			Qlsub=(QFlsub/DENSITY)*pow(10,6);
		    Qrver=(QFlver/DENSITY)*pow(10,6); /*conversion from kg/s to ml/s */
			Qrcar=(QFlcar/DENSITY)*pow(10,6);
			Qrsub=(QFlsub/DENSITY)*pow(10,6);
			
		
			Qao=(QFao/DENSITY)*pow(10,6);
			Qheart=(QFheart/DENSITY)*pow(10,6);
			
			Qlcor1=(QFlcor1/DENSITY)*pow(10,6);
			Qlcor2=(QFlcor1/DENSITY)*pow(10,6);
			Qrcor=(QFlcor1/DENSITY)*pow(10,6);
			
			Qub=Qlver+Qlcar+Qlsub+Qrver+Qrcar+Qrsub;
			Qlb=Qao;
			Qcor=Qlcor1+Qlcor2+Qrcor;

			if (FLAG==1)	/*go to next time-step*/
				for (i=0; i<nSV; i++)
				{
					X[i]=X2[i];
					V[i]=V2[i];
				}
			else		/*decrease time-step*/
				for (i=0; i<nSV; i++)
				{
					X[i]=X1[i];
					V[i]=V1[i];
					
				}
				
		modulo=fmod(T,0.8);
		
		/*Calculation of ev*/
		if (modulo<Ts1)	
		{
		ev=(1-cos((modulo/Ts1)*M_PI))/2;
		}	
		if (modulo>=Ts1&&modulo<Ts2)
		{
		ev=(1+cos((modulo-Ts1)/(Ts2-Ts1)*M_PI))/2;   
		}
		if (modulo>=Ts2)
		{
		ev=0;
		}
	 
		/*Calculation of ea*/
		if (modulo<Tsa1)	
		{
		ea=(1-cos((modulo/Tsa1)*M_PI))/2;
		}	
		if (modulo>=Tsa1&&modulo<Tsa2)
		{
		ea=(1+cos((modulo-Tsa1)/(Tsa2-Tsa1)*M_PI))/2;   
		}
		if (modulo>=Tsa2)
		{
		ea=0;
		}

	/*Equations for SartUB*/
	
	FF[1]=(X[1]-X[2])/R[1];
	DX[1]=(Qub-FF[1])/(C[1]);

	/*Equations for SvnUB*/
	FF[2]=(X[2]-PP[5])/R[2];
	DX[2]=(FF[1]-FF[2])/C[2];

	/*Equations for SartLB*/
	FF[3]=(X[3]-X[4])/R[3];
	DX[3]=(Qlb-FF[3])/(C[3]);
	
	
	/*Equations for SvnLB*/
	FF[4]=(X[4]-PP[5])/R[4];
	DX[4]=(FF[3]-FF[4])/C[4];


	/*Equations for RA*/
	  	  
	 if ((PP[5]-PP[6])>0.0001)
	  FF[5]=sqrt(PP[5]-PP[6])*R[5];
	  else
	  FF[5]=0;
  
	PP[5]=(V[5])*(ea*(EmaxRA-EdRA)+EdRA);
	DV[5]=(FF[4]+FF[2]+FF[12]-FF[5]);
     
	/*Equations for RV*/
	  if ((PP[6]-X[7])>0.01)
	  FF[6]=sqrt(PP[6]-X[7])*R[6];
	  else
	  FF[6]=0;
  
    PP[6]=(V[6])*ev*EmaxRV+(1-ev)*(RV_Pd_alpha*exp(RV_Pd_kappa*(V[6]))+RV_Pd_beta); /*PRESSURE From Elastance*/
	DV[6]=(FF[5]-FF[6]);
	

	/*Equations for pa */
	FF[7]=(X[7]-X[8])/R[7];
	DX[7]=(FF[6]-FF[7])/C[7];
	
	/*Equations for part */
	FF[8]=(X[8]-X[9])/R[8];
	DX[8]=(FF[7]-FF[8])/C[8];
	
	/*Equations for pvn */
	FF[9]=(X[9]-PP[10])/R[9];
	DX[9]=(FF[8]-FF[9])/C[9];

	/*Equations for LA*/
	  if ((PP[10]-PP[11])>0.01)
	  FF[10]=sqrt(PP[10]-PP[11])*R[10];
	  else
	  FF[10]=0;
  
	PP[10]=(V[10])*(ea*(EmaxLA-EdLA)+EdLA);
	DV[10]=(FF[9]-FF[10]);

	/*Equations for LV*/
	  if ((PP[11]-X[0])>0.01)
	  FF[11]=sqrt(PP[11]-X[0])*R[11];
	  else
	  FF[11]=0;
    
    PP[11]=(V[11])*ev*EmaxLV+(1-ev)*(LV_Pd_alpha*exp(LV_Pd_kappa*(V[11]))+LV_Pd_beta); /*PRESSURE From Elastance*/
	DV[11]=(FF[10]-FF[11]);

    /*Equations for AorticSinusCompliace*/
    DX[0]=(FF[11]-Qheart)/0.3;
	
	/*Equations for COR*/
	
	FF[12]=(X[12]-0.75*PP[11])/R[12];
	DX[12]=(Qcor-FF[12])/(C[12]);



/*Implicit Euler for al compartments*/
	X[1]=(X[1]+DT*(Qub/(C[1])+X[2]/(R[1]*(C[1]))))/(1+DT*(1/(R[1]*(C[1]))));
	X[2]=(X[2]+DT*(FF[1]/C[2]+PP[5]/(R[2]*C[2])))/(1+DT*(1/(R[2]*C[2])));
	X[3]=(X[3]+DT*(Qlb/(C[3])+X[4]/(R[3]*(C[3]))))/(1+DT*(1/(R[3]*(C[3]))));
	X[4]=(X[4]+DT*(FF[3]/C[4]+PP[5]/(R[4]*C[4])))/(1+DT*(1/(R[4]*C[4])));
	X[7]=(X[7]+DT*(FF[6]/C[7]+X[8]/(R[7]*C[7])))/(1+DT*(1/(R[7]*C[7])));
	X[8]=(X[8]+DT*(FF[7]/C[8]+X[9]/(R[8]*C[8])))/(1+DT*(1/(R[8]*C[8])));
	X[9]=(X[9]+DT*(FF[8]/C[9]+PP[10]/(R[9]*C[9])))/(1+DT*(1/(R[9]*C[9])));
	
/*Explicit Euler for Heart*/
	V[5]=V[5]+DT*DV[5]; /*RA*/
	V[6]=V[6]+DT*DV[6]; /*RV*/

	V[10]=V[10]+DT*DV[10]; /*%LA*/
	V[11]=V[11]+DT*DV[11]; /*LV*/
	
	X[0]=X[0]+DT*DX[0];

/*Explicit Euler for Coronaries*/
	X[12]=X[12]+DT*DX[12];
	
	
	if (FLAG==1)
	for (i=0; i<nSV; i++)
	{
	X1[i]=X2[i];
	X2[i]=X[i];
	V1[i]=V2[i];
	V2[i]=V[i];

	}
	else
	for (i=0; i<nSV; i++)
	{
	X2[i]=X[i];
	V2[i]=V[i];

	}			
			
			
			
	Pub=(X[1]+Zub*Qub)*133.32; /*conversion from mmHg to Pa Characteristic Impedance is Added!*/
	Plb=(X[3]+Zlb*Qlb)*133.32;
	Pcor=(X[12]+Zcor*Qcor)*133.32;
	Pheart=X[0]*133.32;
	
	
	
	pv=fopen("pv.txt","a");	
	fprintf(pv,"%e %e\n",V[11], PP[11]); 
	fclose(pv);
	

	results=fopen("results.txt","a");
	fprintf(results,"%lf %e %e %e %e %e %e %e %e \n",T,Qheart,Qcor,Qao,(Qlcar+Qlver+Qrcar+Qrver),Pub,Plb,Pcor,Pheart);
	fclose(results);
	

		
			
	}	/*end if(FLAG!=2)*/
		
	
		
		
		
	#endif	/*!RP_NODE*/
	
	host_to_node_real_4(Pub,Plb,Pcor,Pheart);

}	/*end DEFINE_ADJUST*/




DEFINE_PROFILE(inlet_heart, thread, nv)
{
	face_t f=0;
	real QF_par=0.0;
	
	 
	#if !RP_HOST	/*SERIAL or NODE*/
		begin_f_loop(f,thread)
			if(PRINCIPAL_FACE_P(f,thread))
			{	
				QF_par -= (F_FLUX(f,thread));
				
				F_PROFILE(f,thread,nv) = Pheart;
			  
			}
		end_f_loop(f,thread)
	#endif	/*!RP_HOST*/
	
	QFheart=(PRF_GRSUM1(QF_par));
}



DEFINE_PROFILE(outlet_lver, thread, nv)
{
	face_t f=0;
	real QF_par=0.0;
	#if !RP_HOST	/*SERIAL or NODE*/
		begin_f_loop(f,thread)
			if(PRINCIPAL_FACE_P(f,thread))
			{	
				QF_par += (F_FLUX(f,thread));
				F_PROFILE(f,thread,nv) = Pub;
			}
		end_f_loop(f,thread)
	#endif	/*!RP_HOST*/
	
	QFlver=PRF_GRSUM1(QF_par);
}



DEFINE_PROFILE(outlet_lcar, thread, nv)
{
	face_t f=0;
	real QF_par=0.0;
	#if !RP_HOST	/*SERIAL or NODE*/
		begin_f_loop(f,thread)
			if(PRINCIPAL_FACE_P(f,thread))
			{	
				QF_par += (F_FLUX(f,thread));
				F_PROFILE(f,thread,nv) = Pub;
			}
		end_f_loop(f,thread)
	#endif	/*!RP_HOST*/
	
	QFlcar=PRF_GRSUM1(QF_par);
}


DEFINE_PROFILE(outlet_lsub, thread, nv)
{
	face_t f=0;
	real QF_par=0.0;
	#if !RP_HOST	/*SERIAL or NODE*/
		begin_f_loop(f,thread)
			if(PRINCIPAL_FACE_P(f,thread))
			{	
				QF_par += (F_FLUX(f,thread));
				F_PROFILE(f,thread,nv) = Pub;
			}
		end_f_loop(f,thread)
	#endif	/*!RP_HOST*/
	
	QFlsub=PRF_GRSUM1(QF_par);
}



DEFINE_PROFILE(outlet_rver, thread, nv)
{
	face_t f=0;
	real QF_par=0.0;
	#if !RP_HOST	/*SERIAL or NODE*/
		begin_f_loop(f,thread)
			if(PRINCIPAL_FACE_P(f,thread))
			{	
				QF_par += (F_FLUX(f,thread));
				F_PROFILE(f,thread,nv) = Pub;
			}
		end_f_loop(f,thread)
	#endif	/*!RP_HOST*/
	
	QFrver=PRF_GRSUM1(QF_par);
}


DEFINE_PROFILE(outlet_rcar, thread, nv)
{
	face_t f=0;
	real QF_par=0.0;
	#if !RP_HOST	/*SERIAL or NODE*/
		begin_f_loop(f,thread)
			if(PRINCIPAL_FACE_P(f,thread))
			{	
				QF_par += (F_FLUX(f,thread));
				F_PROFILE(f,thread,nv) = Pub;
			}
		end_f_loop(f,thread)
	#endif	/*!RP_HOST*/
	
	QFrcar=PRF_GRSUM1(QF_par);
}


DEFINE_PROFILE(outlet_rsub, thread, nv)
{
	face_t f=0;
	real QF_par=0.0;
	#if !RP_HOST	/*SERIAL or NODE*/
		begin_f_loop(f,thread)
			if(PRINCIPAL_FACE_P(f,thread))
			{	
				QF_par += (F_FLUX(f,thread));
				F_PROFILE(f,thread,nv) = Pub;
			}
		end_f_loop(f,thread)
	#endif	/*!RP_HOST*/
	
	QFrsub=PRF_GRSUM1(QF_par);
}



DEFINE_PROFILE(outlet_aorta, thread, nv)
{
	face_t f=0;
	real QF_par=0.0;
	#if !RP_HOST	/*SERIAL or NODE*/
		begin_f_loop(f,thread)
			if(PRINCIPAL_FACE_P(f,thread))
			{	
				QF_par += (F_FLUX(f,thread));
				F_PROFILE(f,thread,nv) = Plb;
			}
		end_f_loop(f,thread)
	#endif	/*!RP_HOST*/
	
	QFao=PRF_GRSUM1(QF_par);
}


DEFINE_PROFILE(outlet_lcor1, thread, nv)
{
	face_t f=0;
	real QF_par=0.0;
	#if !RP_HOST	/*SERIAL or NODE*/
		begin_f_loop(f,thread)
			if(PRINCIPAL_FACE_P(f,thread))
			{	
				QF_par += (F_FLUX(f,thread));
				F_PROFILE(f,thread,nv) = Pcor;
			}
		end_f_loop(f,thread)
	#endif	/*!RP_HOST*/
	
	QFlcor1=PRF_GRSUM1(QF_par);
}

DEFINE_PROFILE(outlet_lcor2, thread, nv)
{
	face_t f=0;
	real QF_par=0.0;
	#if !RP_HOST	/*SERIAL or NODE*/
		begin_f_loop(f,thread)
			if(PRINCIPAL_FACE_P(f,thread))
			{	
				QF_par += (F_FLUX(f,thread));
				F_PROFILE(f,thread,nv) = Pcor;
			}
		end_f_loop(f,thread)
	#endif	/*!RP_HOST*/
	
	QFlcor2=PRF_GRSUM1(QF_par);
}

DEFINE_PROFILE(outlet_rcor, thread, nv)
{
	face_t f=0;
	real QF_par=0.0;
	#if !RP_HOST	/*SERIAL or NODE*/
		begin_f_loop(f,thread)
			if(PRINCIPAL_FACE_P(f,thread))
			{	
				QF_par += (F_FLUX(f,thread));
				F_PROFILE(f,thread,nv) = Pcor;
			}
		end_f_loop(f,thread)
	#endif	/*!RP_HOST*/
	
	QFrcor=PRF_GRSUM1(QF_par);
}


