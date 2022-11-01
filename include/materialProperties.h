template<int dim>
class properties{
 public:
  properties();
  ~properties();
  void evaluate();
  //private:
  
  void elasticModulus();
  void rotatedElasticModulus();
  void assignGrainAngle();
  void testPrint();
  void crystalRotation();
  std::vector<double> grainAngle;
  std::vector<Table<4, double> > A_phi;
  Table<4, double>ElasticModulus;
  std::vector<FullMatrix<double> > rotationMatrices;
  
};

template<int dim>
properties<dim>::properties():ElasticModulus(dim, dim, dim, dim){

  if(dim ==2){
    double nu12=0.253; double nu21=0.253; double mu=500.0;
    ElasticModulus[0][0][0][0]=alpha1/(1-nu12*nu21);  ElasticModulus[0][0][1][1]=alpha1*nu21/(1-nu12*nu21);
    ElasticModulus[1][1][0][0]=beta1*nu12/(1-nu12*nu21);  ElasticModulus[1][1][1][1]=beta1/(1-nu12*nu21);
    ElasticModulus[0][1][0][1]=mu/2.;   ElasticModulus[0][1][1][0]=mu/2.;
    ElasticModulus[1][0][0][1]=mu/2.;   ElasticModulus[1][0][1][0]=mu/2.;
  }

  if(dim==3){
    double Ex=500.0; double Ey=250.0; double Ez=100.0;
    double nuxy=0.3; double nuyz=0.3; double nuzx=0.3;
    double nuyx=(Ey/Ex)*nuxy; double nuzy=(Ez/Ey)*nuyz; double nuxz=(Ex/Ez)*nuzx;
    double Delta=(1.0-nuxy*nuyx-nuyz*nuzy-nuzx*nuxz-2.0*nuxy*nuyz*nuzx)/(Ex*Ey*Ez);
    double muxy=100.0; double muyz=100.0; double muzx=100.0;
    ElasticModulus[0][0][0][0]=(1.0-nuyz*nuzy)/(Ey*Ez*Delta);
    ElasticModulus[0][0][1][1]=(nuyx+nuzx*nuyz)/(Ey*Ez*Delta);
    ElasticModulus[0][0][2][2]=(nuzx+nuyx*nuzy)/(Ey*Ez*Delta);
    ElasticModulus[1][1][0][0]=(nuxy+nuxz*nuzy)/(Ez*Ex*Delta);
    ElasticModulus[1][1][1][1]=(1.0-nuzx*nuxz)/(Ez*Ex*Delta);
    ElasticModulus[1][1][2][2]=(nuzy+nuzx*nuxy)/(Ez*Ex*Delta);
    ElasticModulus[2][2][0][0]=(nuxz+nuxy*nuyz)/(Ex*Ey*Delta);
    ElasticModulus[2][2][1][1]=(nuyz+nuxz*nuyx)/(Ex*Ey*Delta);
    ElasticModulus[2][2][2][2]=(1.0-nuxy*nuyx)/(Ez*Ex*Delta);
    ElasticModulus[0][1][0][1]=muxy; ElasticModulus[0][1][1][0]=muxy;
    ElasticModulus[0][2][0][2]=muzx; ElasticModulus[0][2][2][0]=muzx;
    ElasticModulus[1][0][0][1]=muxy; ElasticModulus[1][0][1][0]=muxy;
    ElasticModulus[1][2][1][2]=muyz; ElasticModulus[1][2][2][1]=muyz;
    ElasticModulus[2][0][0][2]=muzx; ElasticModulus[2][0][2][0]=muzx;
    ElasticModulus[2][1][1][2]=muyz; ElasticModulus[2][1][2][1]=muyz;
  }
  
}

template<int dim>
properties<dim>::~properties(){}


template<int dim>
void properties<dim>::assignGrainAngle(){
  std::srand(5.5);
  for(unsigned int i=0;i<n_diff_grains;i++){
    double angle= std::rand()%90;
    grainAngle.push_back(angle);
  }

}

template<int dim>
void properties<dim>:: testPrint(){
  //std::cout<<"grain angles:\n";
  /*for(unsigned int i=0;i<n_diff_grains;i++){
    std::cout<<grainAngle[i]<<" ";
  }*/
  /*std::cout<<"\n\nElasticModulus\n\n";
  for(unsigned int i=0;i<dim;i++){
      for(unsigned int j=0;j<dim;j++){
        for(unsigned int k=0;k<dim;k++){
          for(unsigned int l=0;l<dim;l++){
            std::cout<<ElasticModulus[i][j][k][l]<<"  ";
          }std::cout<<"\t\t";
        }
      }
    }
  
  std::cout<<"\n\nA_phi[1]:\n";
  for(unsigned int i=0;i<dim;i++){
    for(unsigned int j=0;j<dim;j++){
      for(unsigned int k=0;k<dim;k++){
	for(unsigned int l=0;l<dim;l++){
	  std::cout<<A_phi[1][i][j][k][l]<<"  ";
	}std::cout<<"\t\t";
      }
    }
  }*/
  //exit(-1);//
}

template<int dim>
void properties<dim>::evaluate(){
  assignGrainAngle();
  //testPrint();
  crystalRotation();
  rotatedElasticModulus();
  testPrint();
}

template<int dim>
void properties<dim>::crystalRotation(){
  FullMatrix<double> rotation(dim, dim);
  for(unsigned int i=0;i<n_diff_grains;i++){
    double var=grainAngle[i]*PI/180.;
    rotation(0,0)=cos(var); rotation(0,1)=-sin(var);
    rotation(1,0)=sin(var);  rotation(1,1)=cos(var);
    if(dim==3){
      rotation(2,2)=1.0;
    }
    rotationMatrices.push_back(rotation);
  }
}


template<int dim>
void properties<dim>:: rotatedElasticModulus(){

  for(unsigned int N=0;N<n_diff_grains;N++){
    Table<4, double>temp(dim, dim, dim, dim);
    FullMatrix<double> rotation(dim, dim);
    rotation=rotationMatrices[N];
    for(unsigned int i=0;i<dim;i++)for(unsigned int j=0;j<dim;j++)for(unsigned int k=0;k<dim;k++)for(unsigned int l=0;l<dim;l++)temp[i][j][k][l]=0.;
    for(unsigned int I=0;I<dim;I++){
      for(unsigned int i=0;i<dim;i++){
        for(unsigned int J=0;J<dim;J++){
          for(unsigned int j=0;j<dim;j++){
            for(unsigned int K=0;K<dim;K++){
              for(unsigned int k=0;k<dim;k++){
                for(unsigned int L=0;L<dim;L++){
                  for(unsigned int l=0;l<dim;l++){
		    temp[i][j][k][l]+=rotation[i][I]*rotation[j][J]*rotation[k][K]*rotation[l][L]*ElasticModulus[I][J][K][L];
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    A_phi.push_back(temp);

  }//n_diff_grain loop ends
  
}
