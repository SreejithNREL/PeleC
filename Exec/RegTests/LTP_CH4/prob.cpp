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

  amrex::Real mw_N2 = 28.0; // N2 molecular weight [g/mol]
  amrex::Real Ru = pele::physics::Constants::RU/mw_N2; // [erg/g/K]
  amrex::Real theta_v_s = 3371.; // [K]
  amrex::Real theta_dof = 5.;
  
  //compute trans-rotational temperature for N2
  amrex::Real Ttr = e/theta_dof/0.5/Ru;
  //get Evib in equilibrium with Ttr
  amrex::Real Evib = Ru*theta_v_s / (std::exp(theta_v_s/Ttr)-1);

  (PeleC::h_prob_parm_device->fuel_state)[URHO]  = rho;
  (PeleC::h_prob_parm_device->fuel_state)[UMX]   = rho * vt;
  (PeleC::h_prob_parm_device->fuel_state)[UMY]   = 0.0;
  (PeleC::h_prob_parm_device->fuel_state)[UMZ]   = 0.0;
  (PeleC::h_prob_parm_device->fuel_state)[UEINT] = rho * e;
  (PeleC::h_prob_parm_device->fuel_state)[UEDEN] = rho * (e + ek);
  (PeleC::h_prob_parm_device->fuel_state)[UTEMP] = T;
  (PeleC::h_prob_parm_device->fuel_state)[UFA]   = rho * Evib;
  for (int n = 0; n < NUM_SPECIES; n++){
    (PeleC::h_prob_parm_device->fuel_state)[UFS + n] = rho * massfrac[n];
  }

}

void
init_composition(amrex::Real molefrac[NUM_SPECIES])
{
  amrex::Real a;
  a = 2.;
  for (int n = 0; n < NUM_SPECIES; n++)
    molefrac[n] = 0.0;

  molefrac[O2_ID]  = 1.0 / (1.0 + PeleC::h_prob_parm_device->phi_in / a + 0.79 / 0.21);
  molefrac[CH4_ID] = PeleC::h_prob_parm_device->phi_in * molefrac[O2_ID] / a;
  molefrac[N2_ID]  = 1.0 - molefrac[CH4_ID] - molefrac[O2_ID];


 // printf("molefrac[CH4_ID] = %f \n", molefrac[CH4_ID]);
 // printf("molefrac[O2_ID] = %f \n", molefrac[O2_ID]);
 // printf("molefrac[N2_ID] = %f \n", molefrac[N2_ID]);
 // amrex::Abort();
  
  // for (int n = 0; n < NUM_SPECIES; n++)
  //   molefrac[n] = 1.0;

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
      pp.query("dilution"       , PeleC::h_prob_parm_device->dilution);
      pp.query("T_in"           , PeleC::h_prob_parm_device->T_in);
      pp.query("vn_in"          , PeleC::h_prob_parm_device->vn_in);
      pp.query("turbulence"     , PeleC::h_prob_parm_device->turbulence);
      pp.query("restart"        , PeleC::h_prob_parm_device->restart);

      pp.query("base_power"     , PeleC::h_prob_parm_device->base_power);
      pp.query("npulses"        , PeleC::h_prob_parm_device->npulses);
      pp.query("pulse_duration" , PeleC::h_prob_parm_device->pulse_duration);
      pp.query("dwell"          , PeleC::h_prob_parm_device->dwell);

      init_bc();

      amrex::Print() << "Base power used for each ignition pulse: " << PeleC::h_prob_parm_device->base_power << std::endl;

//       if (pp.countval("turb_file") > 0) {
// #if AMREX_SPACEDIM==2
//     amrex::Abort("Turbulence inflow unsupported in 2D.");
// #endif

//         PeleC::prob_parm_host->do_turb = true;

//         std::string turb_file = "";
//         pp.query("turb_file", turb_file);
//         amrex::Real turb_scale_loc = 1.0;
//         pp.query("turb_scale_loc", turb_scale_loc);
//         amrex::Real turb_scale_vel = 1.0;
//         pp.query("turb_scale_vel", turb_scale_vel);
//         pp.query("meanFlowDir", PeleC::h_prob_parm_device->meanFlowDir);

//         // Hold nose here - required because of dynamically allocated data in tp
//         AMREX_ASSERT_WITH_MESSAGE(
//         PeleC::h_prob_parm_device->tp.tph == nullptr,
//           "Can only be one TurbParmHost");
//         PeleC::h_prob_parm_device->tp.tph = new TurbParmHost;

//         amrex::Vector<amrex::Real> turb_center = {
//         {0.5 * (probhi[0] + problo[0]), 0.5 * (probhi[1] + problo[1])}};
//         pp.queryarr("turb_center", turb_center);
//         AMREX_ASSERT_WITH_MESSAGE(turb_center.size() == 2, "turb_center must have two elements");
//         for (int n = 0; n < turb_center.size(); ++n) {
//           turb_center[n] *= turb_scale_loc;
//         }

//         amrex::Real restart_time;
//         pp.query("restart_time", restart_time);
//         int turb_nplane = 0;
//         pp.query("turb_nplane", turb_nplane);
//         AMREX_ASSERT(TurbParm->nplane > 0);

//         amrex::Real turb_conv_vel = 0;
//         pp.query("turb_conv_vel", turb_conv_vel);
//         AMREX_ASSERT(TurbParm->turb_conv_vel > 0);

//         init_turbinflow(
//         turb_file, turb_scale_loc, turb_scale_vel, turb_center, turb_conv_vel,restart_time,
//         turb_nplane, PeleC::h_prob_parm_device->tp);
//       }
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

// void PeleC::problem_post_init(){}
void PeleC::problem_post_init(
 int i,
 int j,
 int k,
 amrex::Array4<amrex::Real> const& state,
 amrex::GeometryData const& geomdata){}
//{
//
//  amrex::Real rho = h_prob_parm_device->fuel_state[URHO];
//  amrex::Real velx = state(i, j, k, UMX); //state() at this point is storing velocity, not momentum
//  amrex::Real vely = state(i, j, k, UMY); //state() at this point is storing velocity, not momentum
//  amrex::Real velz = state(i, j, k, UMZ); //state() at this point is storing velocity, not momentum
//
//  amrex::Real sum = 0.0;
//  for (int n = 0; n < NUM_SPECIES; n++){
//    state(i, j, k, UFS + n) = h_prob_parm_device->fuel_state[UFS + n];
//  }
//
//  amrex::Real e = h_prob_parm_device->fuel_state[UEINT];
//
//  state(i, j, k, URHO)  = rho;
//  state(i, j, k, UMX)   = rho * velx;
//  state(i, j, k, UMY)   = rho * vely;
//  state(i, j, k, UMZ)   = rho * velz;
//  state(i, j, k, UTEMP) = h_prob_parm_device->fuel_state[UTEMP];
//  state(i, j, k, UEINT) = e;
//  state(i, j, k, UEDEN) = e + rho * (0.5 * (velx * velx + vely * vely + velz * velz)); 
//
//
//}

void PeleC::problem_post_restart(){}

