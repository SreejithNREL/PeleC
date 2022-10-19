#include "prob.H"
#include <turbinflow.H>

void
init_bc()
{
  amrex::Real vt, ek, a, yl, yr, sumY, T, rho, e;
  amrex::Real molefrac[NUM_SPECIES], massfrac[NUM_SPECIES];

  init_composition(molefrac);


  T = PeleC::h_prob_parm_device->T_in;

  double p = PeleC::h_prob_parm_device->pamb;
  auto eos = pele::physics::PhysicsType::eos();

  eos.X2Y(molefrac, massfrac);
  eos.PYT2RE(p, massfrac, T, rho, e);

  vt = PeleC::h_prob_parm_device->vn_in;
  ek = 0.5 * (vt * vt);

  (PeleC::h_prob_parm_device->fuel_state)[URHO]  = rho;
  (PeleC::h_prob_parm_device->fuel_state)[UMX]   = rho * vt;
  (PeleC::h_prob_parm_device->fuel_state)[UMY]   = 0.0;
  (PeleC::h_prob_parm_device->fuel_state)[UMZ]   = 0.0;
  (PeleC::h_prob_parm_device->fuel_state)[UEINT] = rho * e;
  (PeleC::h_prob_parm_device->fuel_state)[UEDEN] = rho * (e + ek);
  (PeleC::h_prob_parm_device->fuel_state)[UTEMP] = T;
  for (int n = 0; n < NUM_SPECIES; n++)
    (PeleC::h_prob_parm_device->fuel_state)[UFS + n] = rho * massfrac[n];
    // (PeleC::h_prob_parm_device->fuel_state)[UFS + n] = massfrac[n];

}

void
init_composition(amrex::Real molefrac[NUM_SPECIES])
{
  amrex::Real vt, ek, a, yl, yr, sumY, T, rho, e;

   for (int n = 0; n < NUM_SPECIES; n++)
     molefrac[n] = 0.0;
  // molefrac[O2_ID] = 0.21;
  // molefrac[N2_ID] = 0.79;

   // number of carbon and hydrogen atoms in each fuel component
   //order: IC8H18; NC7H16_ID; C6H5CH3; C2H5OH
   amrex::Real C_atoms[4] = {8.0 ,7.0 ,7.0,2.0}; 
   amrex::Real H_atoms[4] = {18.0,16.0,4.0,6.0}; 
   amrex::Real O_atoms[4] = {0.0 ,0.0 ,0.0,1.0}; 

   amrex::Real volfrac_IC8H18  = 0.284;
   amrex::Real volfrac_NC7H16  = 0.172;
   amrex::Real volfrac_C6H5CH3 = 0.339;
   amrex::Real volfrac_C2H5OH  = 0.205;
   
   amrex::Real volfrac_CO2 =  C_atoms[0]*volfrac_IC8H18 + C_atoms[1]*volfrac_NC7H16 + C_atoms[2]*volfrac_C6H5CH3 + C_atoms[3]*volfrac_C2H5OH;
   amrex::Real volfrac_H2O = (H_atoms[0]*volfrac_IC8H18 + H_atoms[1]*volfrac_NC7H16 + H_atoms[2]*volfrac_C6H5CH3 + H_atoms[3]*volfrac_C2H5OH)/2.0;
   
   amrex::Real volfrac_O2  = (2.0*volfrac_CO2+volfrac_H2O-(O_atoms[0]*volfrac_IC8H18  + 
                                                           O_atoms[1]*volfrac_NC7H16  + 
                                                           O_atoms[2]*volfrac_C6H5CH3 + 
                                                           O_atoms[3]*volfrac_C2H5OH))/2.0/PeleC::h_prob_parm_device->phi_in;
   amrex::Real volfrac_N2  = volfrac_O2*(79./21.);

   amrex::Real sum_oxi  = volfrac_IC8H18 + volfrac_NC7H16 + volfrac_C6H5CH3 + volfrac_C2H5OH + volfrac_O2  + volfrac_N2;
   amrex::Real sum_prod = volfrac_CO2    + volfrac_H2O    + volfrac_N2;

   amrex::Real X_IC8H18  = volfrac_IC8H18   / sum_oxi;
   amrex::Real X_NC7H16  = volfrac_NC7H16   / sum_oxi;
   amrex::Real X_C6H5CH3 = volfrac_C6H5CH3  / sum_oxi;
   amrex::Real X_C2H5OH  = volfrac_C2H5OH   / sum_oxi;
   amrex::Real X_O2      = volfrac_O2       / sum_oxi;
   amrex::Real X_N2      = volfrac_N2       / sum_oxi;

   amrex::Real X_co2    = volfrac_CO2 / sum_prod;
   amrex::Real X_h2o    = volfrac_H2O / sum_prod;

   for (int n = 0; n < NUM_SPECIES; n++)
     molefrac[n] = 0.0;

   molefrac[IC8_ID]     = X_IC8H18  * (1.0-PeleC::h_prob_parm_device->egr);
   molefrac[NC7H16_ID]  = X_NC7H16  * (1.0-PeleC::h_prob_parm_device->egr);
   molefrac[C6H5CH3_ID] = X_C6H5CH3 * (1.0-PeleC::h_prob_parm_device->egr);
   molefrac[C2H5OH_ID]  = X_C2H5OH  * (1.0-PeleC::h_prob_parm_device->egr);
   molefrac[O2_ID]      = X_O2      * (1.0-PeleC::h_prob_parm_device->egr);
   molefrac[N2_ID]      = X_N2      * (1.0-PeleC::h_prob_parm_device->egr) + (PeleC::h_prob_parm_device->egr*volfrac_N2 / sum_prod);
   molefrac[CO2_ID]     = X_co2*PeleC::h_prob_parm_device->egr;
   molefrac[H2O_ID]     = X_h2o*PeleC::h_prob_parm_device->egr;
}

void
pc_prob_close()
{

}

extern "C" 
{
  void amrex_probinit(
              const int* init,
              const int* name,
              const int* namelen,
              const amrex_real* problo,
              const amrex_real* probhi)
    {

      amrex::ParmParse pp("prob");
      pp.query("pamb"           , PeleC::h_prob_parm_device->pamb);
      pp.query("phi_in"         , PeleC::h_prob_parm_device->phi_in);
      pp.query("egr"            , PeleC::h_prob_parm_device->egr);
      pp.query("T_in"           , PeleC::h_prob_parm_device->T_in);
      pp.query("vn_in"          , PeleC::h_prob_parm_device->vn_in);
      pp.query("turbulence"     , PeleC::h_prob_parm_device->turbulence);
      pp.query("init_kernel"    , PeleC::h_prob_parm_device->init_kernel);
      pp.query("iname"          , PeleC::prob_parm_host->iname);
      pp.query("time_init_turb" , PeleC::h_prob_parm_device->time_init_turb); 
      pp.query("restart"        , PeleC::h_prob_parm_device->restart);

      PeleC::h_prob_parm_device->L[0] = probhi[0] - problo[0];
      PeleC::h_prob_parm_device->L[1] = probhi[1] - problo[1];
      PeleC::h_prob_parm_device->L[2] = probhi[2] - problo[2];


      if (PeleC::h_prob_parm_device->restart) {
        amrex::Print() << "Skipping input file reading and assuming restart."
                       << std::endl;
      } 
      else if(PeleC::h_prob_parm_device->init_kernel == false){
        //initialize premixed state
        init_bc();
      }
      else {
        //get array sizes
        int z_coord = 0;
        get_inputs(PeleC::prob_parm_host->iname,z_coord);
        
        PeleC::h_prob_parm_device->nx = PeleC::prob_parm_host->nx_tmp;
        PeleC::h_prob_parm_device->ny = PeleC::prob_parm_host->ny_tmp;
        PeleC::h_prob_parm_device->nz = PeleC::prob_parm_host->nz_tmp;
        PeleC::h_prob_parm_device->nscal = PeleC::prob_parm_host->nscal_tmp;

        amrex::Vector<double> data_tmp(PeleC::h_prob_parm_device->nx * PeleC::h_prob_parm_device->ny * PeleC::h_prob_parm_device->nscal);
        amrex::Vector<double> gridx_input(PeleC::h_prob_parm_device->nx);
        amrex::Vector<double> gridy_input(PeleC::h_prob_parm_device->ny);
        amrex::Vector<double> gridz_input(PeleC::h_prob_parm_device->nz);

        // //get data from file
        read_input_file(PeleC::prob_parm_host->iname,
                        z_coord,
                        PeleC::h_prob_parm_device->nx,
                        PeleC::h_prob_parm_device->ny,
                        PeleC::h_prob_parm_device->nz,
                        PeleC::h_prob_parm_device->nscal,
                        gridx_input, gridy_input, gridz_input, data_tmp);

        //copy data to the working arrays
        PeleC::prob_parm_host->v_xinput.resize(PeleC::h_prob_parm_device->nx);
        for (int i = 0; i < PeleC::h_prob_parm_device->nx; i++) {
          PeleC::prob_parm_host->v_xinput[i] = gridx_input[i];
        }

        PeleC::prob_parm_host->v_yinput.resize(PeleC::h_prob_parm_device->ny);
        for (int i = 0; i < PeleC::h_prob_parm_device->ny; i++) {
          PeleC::prob_parm_host->v_yinput[i] = gridy_input[i];
        }

        PeleC::prob_parm_host->data.resize(PeleC::h_prob_parm_device->nx*PeleC::h_prob_parm_device->ny*PeleC::h_prob_parm_device->nscal);
        for (int i = 0; i < PeleC::h_prob_parm_device->nx*PeleC::h_prob_parm_device->ny*PeleC::h_prob_parm_device->nscal; i++) {
          PeleC::prob_parm_host->data[i] = data_tmp[i];
          // printf("data_tmp[i] %f \n", data_tmp[i]);
        }


        // Get the array differences.
        PeleC::prob_parm_host->v_xdiff.resize(PeleC::h_prob_parm_device->nx);
        std::adjacent_difference(
          PeleC::prob_parm_host->v_xinput.begin(), PeleC::prob_parm_host->v_xinput.end(),
          PeleC::prob_parm_host->v_xdiff.begin());
        PeleC::prob_parm_host->v_xdiff[0] = PeleC::prob_parm_host->v_xdiff[1];

        PeleC::prob_parm_host->v_ydiff.resize(PeleC::h_prob_parm_device->ny);
        std::adjacent_difference(
          PeleC::prob_parm_host->v_yinput.begin(), PeleC::prob_parm_host->v_yinput.end(),
          PeleC::prob_parm_host->v_ydiff.begin());
        PeleC::prob_parm_host->v_ydiff[0] = PeleC::prob_parm_host->v_ydiff[1];

        // Get pointer to the data
        PeleC::prob_parm_host->xarray.resize(PeleC::prob_parm_host->v_xinput.size());
        PeleC::prob_parm_host->xdiff.resize(PeleC::prob_parm_host->v_xdiff.size());
        PeleC::prob_parm_host->yarray.resize(PeleC::prob_parm_host->v_yinput.size());
        PeleC::prob_parm_host->ydiff.resize(PeleC::prob_parm_host->v_ydiff.size());
        PeleC::prob_parm_host->data_ic.resize(PeleC::prob_parm_host->data.size());

        amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, PeleC::prob_parm_host->v_xinput.begin(),
         PeleC::prob_parm_host->v_xinput.end(),
         PeleC::prob_parm_host->xarray.begin());
            
        amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, PeleC::prob_parm_host->v_xdiff.begin(),
         PeleC::prob_parm_host->v_xdiff.end(),
         PeleC::prob_parm_host->xdiff.begin());
            
        amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, PeleC::prob_parm_host->v_yinput.begin(),
         PeleC::prob_parm_host->v_yinput.end(),
         PeleC::prob_parm_host->yarray.begin());
            
        amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, PeleC::prob_parm_host->v_ydiff.begin(),
         PeleC::prob_parm_host->v_ydiff.end(),
         PeleC::prob_parm_host->ydiff.begin());
            
        amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, PeleC::prob_parm_host->data.begin(),
         PeleC::prob_parm_host->data.end(),
         PeleC::prob_parm_host->data_ic.begin());

        PeleC::h_prob_parm_device->d_xarray  = PeleC::prob_parm_host->xarray.data();
        PeleC::h_prob_parm_device->d_xdiff   = PeleC::prob_parm_host->xdiff.data();
        PeleC::h_prob_parm_device->d_yarray  = PeleC::prob_parm_host->yarray.data();
        PeleC::h_prob_parm_device->d_ydiff   = PeleC::prob_parm_host->ydiff.data();
        PeleC::h_prob_parm_device->d_data_ic = PeleC::prob_parm_host->data_ic.data();


        // Dimensions of the input box.
        PeleC::h_prob_parm_device->Lxinput =
          PeleC::prob_parm_host->v_xinput[PeleC::h_prob_parm_device->nx - 1] + 0.5 * PeleC::prob_parm_host->v_xdiff[PeleC::h_prob_parm_device->nx - 1];
        PeleC::h_prob_parm_device->Lyinput =
          PeleC::prob_parm_host->v_yinput[PeleC::h_prob_parm_device->ny - 1] + 0.5 * PeleC::prob_parm_host->v_ydiff[PeleC::h_prob_parm_device->ny - 1];

        //initialize premixed state
        init_bc();
      }
    }
}


void
read_input_file(
  const std::string iname,
  int z_coord,
  int nx,
  int ny,
  int nz,
  int nscal,
  amrex::Vector<amrex::Real>& gridx_input,
  amrex::Vector<amrex::Real>& gridy_input,
  amrex::Vector<amrex::Real>& gridz_input,
  amrex::Vector<amrex::Real>& data)
{
  std::ifstream infile(iname, std::ios::in);
  const std::string memfile = read_file(infile);
  if (not infile.is_open()) {
    amrex::Abort("Unable to open input file " + iname);
  }
  infile.close();
  int i,j,k,n;
  
  std::istringstream iss(memfile);

  // Read the file
  int nlines = 0;
  std::string firstline, line;
  if(z_coord){
    for (i=0;i<nscal+7;i++){
      std::getline(iss, firstline); // skip header
    }  
  }
  else{
    for (i=0;i<nscal+5;i++){
      std::getline(iss, firstline); // skip header
    }    
  }


  // Read grid in X
  int cnt = 0;
  for (i=0;i<nx;i++){
    std::getline(iss, line);
    std::istringstream sinput(line);
    sinput >> gridx_input[cnt];

    if(z_coord == 0){
      gridx_input[cnt] = gridx_input[cnt]*100.0; //convert to [cm]
    }
    cnt++;
  }

  // Read the data from the file
  cnt = 0;
  for (j=0;j<ny;j++){
    std::getline(iss, line);
    std::istringstream sinput(line);
    sinput >> gridy_input[cnt];
    if(z_coord == 0){
      gridy_input[cnt] = gridy_input[cnt]*100.0; //convert to [cm]
    }
    cnt++;
  }
  
  if(z_coord){
    // Read the data from the file
    cnt = 0;
    for (k=0;k<nz;k++){
      std::getline(iss, line);
      std::istringstream sinput(line);
      sinput >> gridz_input[cnt];
      // gridz_input[cnt] = gridz_input[cnt]; //convert to [cm]
      cnt++;
    }
  }

  // Read the data from the file
  cnt = 0;
  //while (std::getline(iss, line)) {
  for (n=0;n<nscal;n++){
    for (k=0;k<nz;k++){
      for (j=0;j<ny;j++){
        for (i=0;i<nx;i++){
          std::getline(iss, line);
          std::istringstream sinput(line);
          sinput >> data[cnt]; //read only the specified number of points in x
          cnt++;
        }
      }
    }
  }
  //if(cnt != PeleC::h_prob_parm_device->nscal*PeleC::h_prob_parm_device->nx*PeleC::h_prob_parm_device->ny){
  //  printf("Number of lines = %i. Specified in the header = %i \n", cnt,PeleC::h_prob_parm_device->nscal*PeleC::h_prob_parm_device->nx*PeleC::h_prob_parm_device->ny);
  //  amrex::Abort("Number of lines in file differ from specified in header");
  //}
};

void
get_inputs(
  const std::string& iname,
  int z_coord)
{
  std::ifstream infile(iname, std::ios::in);
  const std::string memfile = read_file(infile); //The problem is coming from here
  if (not infile.is_open()) {
    amrex::Abort("Unable to open input file " + iname);
  }
  infile.close();
  std::istringstream iss(memfile);

  // Read the file
  std::string line;
  
  std::getline(iss, line);
  std::istringstream ninput(line);
  ninput >> PeleC::prob_parm_host->nscal_tmp;

  std::getline(iss, line);
  std::istringstream xinput(line);
  xinput >> PeleC::prob_parm_host->nx_tmp;
  
  std::getline(iss, line);
  std::istringstream yinput(line);
  yinput >> PeleC::prob_parm_host->ny_tmp;

  amrex::Print() << "Reading data from: " <<  iname.c_str() << std::endl;
  if(z_coord){
    std::getline(iss, line);
    std::istringstream zinput(line);
    zinput >> PeleC::prob_parm_host->nz_tmp;
    PeleC::prob_parm_host->nscal_tmp = PeleC::prob_parm_host->nscal_tmp - 3; //remove x and y and z from count
    amrex::Print() << PeleC::prob_parm_host->nscal_tmp << " scalars found with nx =  " << PeleC::prob_parm_host->nx_tmp 
                          << " and ny = " << PeleC::prob_parm_host->ny_tmp << " and nz = " << PeleC::prob_parm_host->nz_tmp << " points \n" <<std::endl;
  }
  else{
    PeleC::prob_parm_host->nscal_tmp = PeleC::prob_parm_host->nscal_tmp - 2; //remove x and y from count
    PeleC::prob_parm_host->nz_tmp = 1;
    amrex::Print() << PeleC::prob_parm_host->nscal_tmp << " scalars found with nx =  " << PeleC::prob_parm_host->nx_tmp 
                          << " and ny = " << PeleC::prob_parm_host->ny_tmp << " points \n" << std::endl;
  }

};


void PeleC::problem_post_timestep()
{
}

// void PeleC::problem_post_init()
// {
// }

void PeleC::problem_post_restart(){} 

// void PeleC::problem_post_init(){};  

void PeleC::problem_post_init(//){};  
 int i,
 int j,
 int k,
 amrex::Array4<amrex::Real> const& state,
 amrex::GeometryData const& geomdata)
{
 const amrex::Real* prob_lo = geomdata.ProbLo();
 const amrex::Real* prob_hi = geomdata.ProbHi();
 const amrex::Real* dx = geomdata.CellSize();
 const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
 const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
 const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
 amrex::Real u[3] = {0.0};
 amrex::Real molefrac[NUM_SPECIES] = {0.0};
 amrex::Real massfrac[NUM_SPECIES] = {0.0};
 amrex::Real rho, energy, T, pres, e;

 amrex::Real kernel_location[3];
 const amrex::Real kernel_diameter = 100.0e-04;
 const amrex::Real kernel_y_loc    = 1.92;
 const amrex::Real kernel_height   = 1.01e-01;
 double r;

 kernel_location[0] = 0.0;
 kernel_location[1] = kernel_y_loc-kernel_height/2.0;
 kernel_location[2] = 0.0;

 if(PeleC::h_prob_parm_device->init_kernel){
   
   r = sqrt(pow((kernel_location[0]-x),2)+pow((kernel_location[2]-z),2));
   if(y > h_prob_parm_device->d_yarray[0] and y < h_prob_parm_device->d_yarray[PeleC::h_prob_parm_device->ny-1] 
     and r < PeleC::h_prob_parm_device->d_xarray[PeleC::h_prob_parm_device->nx-1]){

     // Fill in the velocities and interpolated quantities.
     amrex::Real u[3] = {0.0};
     amrex::Real uinterp[h_prob_parm_device->nscal] = {0.0};

     // Interpolation factors
     amrex::Real mod[3] = {0.0};
     int idx   = 0;
     int idxp1 = 0;
     int idy   = 0;
     int idyp1 = 0;
     amrex::Real slp[3] = {0.0};
     
     mod[0] = std::fmod(r, h_prob_parm_device->Lxinput);
     locate(h_prob_parm_device->d_xarray, h_prob_parm_device->nx, mod[0], idx);
     idxp1 = (idx + 1) % h_prob_parm_device->nx;
     //idxp1 = (idx + 1);
     slp[0] = (mod[0] - h_prob_parm_device->d_xarray[idx]) / h_prob_parm_device->d_xdiff[idx];

     mod[1] = std::fmod(y, h_prob_parm_device->Lyinput);
     locate(h_prob_parm_device->d_yarray, h_prob_parm_device->ny, mod[1], idy);
     idyp1 = (idy + 1) % h_prob_parm_device->ny;
     //idyp1 = (idy + 1);
     slp[1] = (mod[1] - h_prob_parm_device->d_yarray[idy]) / h_prob_parm_device->d_ydiff[idy];

     const amrex::Real f0 = (1 - slp[0]) * (1 - slp[1]);
     const amrex::Real f1 = slp[0] * (1 - slp[1]);
     const amrex::Real f2 = (1 - slp[0]) * slp[1];
     const amrex::Real f3 = slp[0] * slp[1];

     // Interpolate data
     for (int iscal=0;iscal<h_prob_parm_device->nscal;iscal++){
       uinterp[iscal] = h_prob_parm_device->d_data_ic[idx  +idy  *h_prob_parm_device->nx+iscal*h_prob_parm_device->nx*h_prob_parm_device->ny]*f0 +
                        h_prob_parm_device->d_data_ic[idxp1+idy  *h_prob_parm_device->nx+iscal*h_prob_parm_device->nx*h_prob_parm_device->ny]*f1 +
                        h_prob_parm_device->d_data_ic[idx  +idyp1*h_prob_parm_device->nx+iscal*h_prob_parm_device->nx*h_prob_parm_device->ny]*f2 +
                        h_prob_parm_device->d_data_ic[idxp1+idyp1*h_prob_parm_device->nx+iscal*h_prob_parm_device->nx*h_prob_parm_device->ny]*f3;
     
     }

     // // temperature, pressure, internal energy and composition
     double p = uinterp[1]*10.0; //Pa to dyn/cm2
     // rho = uinterp[0]/1.0e+3;
     T = uinterp[2];
     double massfrac[NUM_SPECIES]={0.0};

     for (int n = 0; n < NUM_SPECIES; n++)
       massfrac[n] = uinterp[n+6];

     double sum = 0.0;
     for (int n = 0; n < NUM_SPECIES; n++)
       sum = sum + massfrac[n];
     // printf("Sum of all species = % e \n", sum);
     massfrac[N2_ID] = massfrac[N2_ID] + 1.0 - sum; //making sure that sum of ys is 1


     auto eos = pele::physics::PhysicsType::eos();
     // eos.RTY2P(rho, T, massfrac, p);
     eos.PYT2RE(p, massfrac, T, rho, e);

     double beta  = atan(z/x);
     double alpha = atan(x/z);

     //Velocity in each direction
     if(x >= 0.0){
       u[0] = uinterp[5]*cos(beta)*100.0;
     }
     else{
       u[0] = -uinterp[5]*cos(beta)*100.0;
     }

     u[1] = uinterp[4]*100.0;

     if(z >= 0){
       u[2] = uinterp[5]*cos(alpha)*100.0;
     }
     else{
       u[2] = -uinterp[5]*cos(alpha)*100.0;
     }


     // amrex::Real rho_init = state(i, j, k, URHO);
     amrex::Real velx = u[0] + state(i, j, k, UMX); //state() at this point is storing velocity, not momentum
     amrex::Real vely = u[1] + state(i, j, k, UMY); //state() at this point is storing velocity, not momentum
     amrex::Real velz = u[2] + state(i, j, k, UMZ); //state() at this point is storing velocity, not momentum

     for (int n = 0; n < NUM_SPECIES; n++)
       state(i, j, k, UFS + n) = rho * massfrac[n];

     state(i, j, k, URHO)  = rho;
     state(i, j, k, UTEMP) = T;
     state(i, j, k, UMX)   = rho * velx;
     state(i, j, k, UMY)   = rho * vely;
     state(i, j, k, UMZ)   = rho * velz;
     state(i, j, k, UEINT) = rho * e;
     state(i, j, k, UEDEN) = rho * (e + 0.5 * (velx * velx + vely * vely + velz * velz)); 

     amrex::Real spec_sum = 0.0;
     for (int n = 0; n < NUM_SPECIES; n++) {
       spec_sum = spec_sum + state(i, j, k, UFS + n);
     }
     if (
       amrex::Math::abs(state(i, j, k, URHO) - spec_sum) >
       1.e-6 * state(i, j, k, URHO)) {
       printf("In problem_post_restart: Sum of (rho X)_i = %e vs rho = %e at (%i,%i,%i)\n",spec_sum,state(i,j,k,URHO),i,j,k);
     }

   }  
   else{//initialize premixed state outside the kernel region
     // amrex::Real rho  = state(i, j, k, URHO);
     amrex::Real velx = state(i, j, k, UMX); //state() at this point is storing velocity, not momentum
     amrex::Real vely = state(i, j, k, UMY); //state() at this point is storing velocity, not momentum
     amrex::Real velz = state(i, j, k, UMZ); //state() at this point is storing velocity, not momentum
     state(i, j, k,URHO)  = h_prob_parm_device->fuel_state[URHO];
     state(i, j, k,UMX)   = h_prob_parm_device->fuel_state[URHO]*velx;
     state(i, j, k,UMY)   = h_prob_parm_device->fuel_state[URHO]*vely;
     state(i, j, k,UMZ)   = h_prob_parm_device->fuel_state[URHO]*velz;
     state(i, j, k,UEINT) = h_prob_parm_device->fuel_state[UEINT];
     state(i, j, k,UEDEN) = h_prob_parm_device->fuel_state[UEINT] + h_prob_parm_device->fuel_state[URHO] * (0.5 * (velx * velx + vely * vely + velz * velz));
     state(i, j, k,UTEMP) = h_prob_parm_device->fuel_state[UTEMP];
     for (int n = 0; n < NUM_SPECIES; n++)
       state(i, j, k,UFS + n) = h_prob_parm_device->fuel_state[UFS + n];


   }
 }
 else{ // init_kernel = false: initialize premixed state with velocity field coming from plt file
   // amrex::Real rho  = state(i, j, k, URHO)
   amrex::Real velx = state(i, j, k, UMX);
   amrex::Real vely = state(i, j, k, UMY);
   amrex::Real velz = state(i, j, k, UMZ);
   state(i, j, k,URHO)  = h_prob_parm_device->fuel_state[URHO];
   state(i, j, k,UMX)   = h_prob_parm_device->fuel_state[URHO]*velx;
   state(i, j, k,UMY)   = h_prob_parm_device->fuel_state[URHO]*vely;
   state(i, j, k,UMZ)   = h_prob_parm_device->fuel_state[URHO]*velz;
   state(i, j, k,UEINT) = h_prob_parm_device->fuel_state[UEINT];
   state(i, j, k,UEDEN) = h_prob_parm_device->fuel_state[UEINT] + h_prob_parm_device->fuel_state[URHO] * (0.5 * (velx * velx + vely * vely + velz * velz));
   state(i, j, k,UTEMP) = h_prob_parm_device->fuel_state[UTEMP];
   for (int n = 0; n < NUM_SPECIES; n++)
     state(i, j, k,UFS + n) = h_prob_parm_device->fuel_state[UFS + n];
 }

}
