//===================================================================
//    ELS Parameter Class : ELSParameters.cc
//      Author    : T.Shibata
//      Creation  : 2009.11.30
//      Last Update ; 2010.11.10
//      Last Update ; 2011.01.24
//      Last Update ; 2011.02.03
//      Last Update ; 2011.03.23
//      Last Update : 2011.04.12  
//      Last Update : 2011.05.15  
//      Last Update : 2011.12.11
//      Last Update : 2011.12.12
//      Last Update : 2011.12.13
//      Last Update : 2011.12.14 
//      Last Update : 2011.12.25
//      Last Update : 2012.01.07
//      Last Update : 2012.01.11 ( Bug Fixed )
//      Last Update : 2012.05.01
//      Last Update : 2012.09.24
//      Last Update : 2012.10.02
//===================================================================
#include "ELSParameters.hh"
#include "PhysicsParameters.hh"
//-----------------------------------------------
using namespace std;
//------------------------------------------------
ELSParameters::ELSParameters(){}
//------------------------------------------------
ELSParameters::~ELSParameters(){}
//------------------------------------------------
void ELSParameters::els_geometory_initialize(Plist *plist)
{

  els_alignment_setting(plist);

  els_parameter_initialize(plist);

  shielding_parameter_initialize();  

  detectorplane_parameter_initialize();

  dump_ELS_geometory();
}
//------------------------------------------------
void ELSParameters::els_alignment_setting(Plist *plist)
{

  elssite_ground_height = plist->ELSSiteGroundHeight(0);

  // Alignment case-1 : Cu Collimeter Hole Position
  alignment_position_CuCollimeter_hole[0]=0.0;  // x-axis has no mean. 
  alignment_position_CuCollimeter_hole[1]=0.0;  // maximum shift is +-1 mm --> 0 mm no alighmen
  alignment_position_CuCollimeter_hole[2]=plist->CuCollimeter_shiftz();  // 
  
  // Alignment case-2 : Vertical Beam line position  
  alignment_vertical_beamline_shift[0] = plist->Vertical_beamline_shift(0);
  alignment_vertical_beamline_shift[1] = plist->Vertical_beamline_shift(1);

  // Alignment case-3 : BM York & Coil Position
  alignment_position_BMYork_Coil[0] = plist->BMYork_Coil_shift(0);
  alignment_position_BMYork_Coil[1] = plist->BMYork_Coil_shift(1);
  alignment_position_BMYork_Coil[2] = plist->BMYork_Coil_shift(2);

  // Alignment case-4 : Slit Collimeter  
  slit_width=plist->SlitWidth();
  alignment_position_Slit_Collimeter[0]=plist->Slit_Collimeter_shift(0); 
  alignment_position_Slit_Collimeter[1]=plist->Slit_Collimeter_shift(1); 

  // Alignment case-5 : Faraday Cup
  alignment_position_Faradaycup[0]=plist->Faradaycup_shift(0);
  alignment_position_Faradaycup[1]=plist->Faradaycup_shift(1);
  alignment_position_Faradaycup[2]=plist->Faradaycup_shift(2);

  // Alignment case-6 : Cover Box Hole
  alignment_coverbox_boardhole_shift[0]=plist->Coverbox_Hole_shift(0);
  alignment_coverbox_boardhole_shift[1]=plist->Coverbox_Hole_shift(1);

}
//------------------------------------------------
void ELSParameters::els_parameter_initialize(Plist *plist)
{
  //-------------------------------------------------
  // ELS Fence &
  // Concrete pad & 
  //  ELS ( 40ft HighCube ) Container // unit = mm
  // ------------------------------------------------

  // Gometory
  size_of_concretepad[0] = 71.4*ft2mm;
  size_of_concretepad[1] = 14.0*ft2mm;
  size_of_concretepad[2] = 10.0*inch2mm;

  outer_size_of_ELScontainer[0] =  12192.0;   
  outer_size_of_ELScontainer[1] =  2438.4;
  outer_size_of_ELScontainer[2] =  2895.6;
  
  distance_concretepad_ELS[0] = size_of_concretepad[0] 
                                -  ( outer_size_of_ELScontainer[0] + 1000.0 );
  distance_concretepad_ELS[1] = (size_of_concretepad[1]  
                                  - ( outer_size_of_ELScontainer[1] + 158.0*2.0 ))/2.0;
  distance_concretepad_ELS[2] = 0.0;
  /* 1000.0 is measurement value between edge of concrete pad to ELS backward wall.
     158.9 is measurement value of thick of insulation wall */

  thick_of_ELScontainer_wall_forward  = 8.0;
  thick_of_ELScontainer_wall_backward = 8.0;
  thick_of_ELScontainer_wall_side     = 8.0;
  thick_of_ELScontainer_wall_roof     = 8.0;
  thick_of_ELScontainer_wall_bottom   = 155.3;

  inner_size_of_ELScontainer[0] = outer_size_of_ELScontainer[0]
                                   - ( thick_of_ELScontainer_wall_forward +
                                       thick_of_ELScontainer_wall_backward );
                                  // true inner length = 12,032.0
  inner_size_of_ELScontainer[1] = outer_size_of_ELScontainer[1]
                                   - ( thick_of_ELScontainer_wall_side * 2.0 );
                                  // true inner width  = 2352.0
  inner_size_of_ELScontainer[2] = outer_size_of_ELScontainer[2]
                                   - ( thick_of_ELScontainer_wall_roof +
                                       thick_of_ELScontainer_wall_bottom );  
                                  // true inner height is 2698.0
 
   position_bias_ELScontainer_wall_forward  = 142.0;
   position_bias_ELScontainer_wall_backward = 50.0;  
   position_bias_ELScontainer_wall_side     = 43.2;  

    // Position
    position_of_concretepad[0] = size_of_concretepad[0]/2.0-distance_concretepad_ELS[0];
    position_of_concretepad[1] = -size_of_concretepad[1]/2.0+distance_concretepad_ELS[1]+158.0;
    position_of_concretepad[2] = elssite_ground_height;  //0.0;

    position_of_ELS_outer[0] = outer_size_of_ELScontainer[0]/2.0;
    position_of_ELS_outer[1] = -outer_size_of_ELScontainer[1]/2.0;
    position_of_ELS_outer[2] = position_of_concretepad[2]
                               + size_of_concretepad[2]/2.0 
                               + outer_size_of_ELScontainer[2]/2.0;
                                 
    position_of_ELS_inner[0] = thick_of_ELScontainer_wall_forward
                             + inner_size_of_ELScontainer[0]/2.0;
    position_of_ELS_inner[1] = -thick_of_ELScontainer_wall_side
                               - inner_size_of_ELScontainer[1]/2.0;
    position_of_ELS_inner[2] = position_of_concretepad[2]
                               + size_of_concretepad[2]/2.0
                               + thick_of_ELScontainer_wall_bottom 
                               + inner_size_of_ELScontainer[2]/2.0;

    position_of_ELS[0]=position_of_ELS_outer[0];
    position_of_ELS[1]=position_of_ELS_outer[1];
    position_of_ELS[2]=position_of_ELS_outer[2];

    for( int i=0; i<3; i++ ){
         rotation_of_concretepad[i]=0.0;   
         rotation_of_ELS[i]=0.0;
         rotation_of_ELS_outer[i]=0.0;
         rotation_of_ELS_inner[i]=0.0;
    }
                                                                    
    //-------------------------------------------------
    //  Base A & B & C
    //  id=0 : A
    //  id=1 : B
    //  id=2 : C
    //-------------------------------------------------

    // Geometory

    size_of_baseA[0] = 6460.0;
    size_of_baseA[1] = 2298.0;
    size_of_baseA[2] = 309.0;

    size_of_baseB[0] = 5140.0;
    size_of_baseB[1] = 2298.0;
    size_of_baseB[2] = 259.0;

    size_of_baseC[0] = outer_size_of_ELScontainer[0]  
                        - position_bias_ELScontainer_wall_forward  
                         - position_bias_ELScontainer_wall_backward;
    size_of_baseC[1] = outer_size_of_ELScontainer[1]
                        - position_bias_ELScontainer_wall_side*2.0;
    size_of_baseC[2] = 9.0;

    space_bwt_base_wall_forward = 200.0;
    space_bwt_base_wall_backward = 200.0;
    space_bwt_base_wall_side = 62.2;  // 62.2 = 70.2 - 8.0

    // Position
    position_of_baseB[0] = position_bias_ELScontainer_wall_forward
                            + space_bwt_base_wall_forward
                             + size_of_baseB[0]/2.0;
    position_of_baseB[1] = - thick_of_ELScontainer_wall_side
                            - space_bwt_base_wall_side 
                             - size_of_baseB[1]/2.0;
    position_of_baseB[2] = position_of_concretepad[2]
                            + size_of_concretepad[2]/2.0
                             + thick_of_ELScontainer_wall_bottom
                              + size_of_baseC[2]
                               + size_of_baseB[2]/2.0;

    position_of_baseA[0] = position_of_baseB[0]
                            + size_of_baseB[0]/2.0
                              + size_of_baseA[0]/2.0;
    position_of_baseA[1] = position_of_baseB[1];
    position_of_baseA[2] = position_of_concretepad[2]
                            + size_of_concretepad[2]/2.0
                             + thick_of_ELScontainer_wall_bottom
                              + size_of_baseC[2]
                               + size_of_baseA[2]/2.0;
  
   position_of_baseC[0] = position_bias_ELScontainer_wall_forward
                           + size_of_baseC[0]/2.0;
   position_of_baseC[1] = - position_bias_ELScontainer_wall_side
                            - size_of_baseC[1]/2.0;
   position_of_baseC[2] = position_of_concretepad[2]
                            + size_of_concretepad[2]/2.0
                             + thick_of_ELScontainer_wall_bottom
                              + size_of_baseC[2]/2.0;

   for( int i=0; i<3; i++ ){
     rotation_of_baseA[i]=0.0;
     rotation_of_baseB[i]=0.0;
     rotation_of_baseC[i]=0.0;
   }

   for( int i=0; i<3; i++ ){
     size_of_base[0][i]=size_of_baseA[i];
     size_of_base[1][i]=size_of_baseB[i];
     size_of_base[2][i]=size_of_baseC[i];
        
     position_of_base[0][i]=position_of_baseA[i];
     position_of_base[1][i]=position_of_baseB[i];
     position_of_base[2][i]=position_of_baseC[i];

     rotation_of_base[0][i]=rotation_of_baseA[i];
     rotation_of_base[1][i]=rotation_of_baseB[i];
     rotation_of_base[2][i]=rotation_of_baseC[i];
   }

   //-------------------------------------------------
   //  Beam Line Component
   //-------------------------------------------------

   position_of_beam_start[0] = position_of_baseA[0] - size_of_baseA[0]/2.0;
   position_of_beam_start[1] = position_of_baseA[1] - size_of_baseA[1]/2.0
                               + 650.0;
   position_of_beam_start[2] = position_of_baseA[2] + size_of_baseA[2]/2.0
                               + 700.0;

   // Frange ICF203(1)
   //RMin_ICF203 = 152.0/2.0;
   RMin_ICF203 = 20.0/2.0;
   RMax_ICF203 = 203.0/2.0;
   L_ICF203    = 22.0;

   // Frange ICF070(1)
   RMin_ICF070_1 = 38.3/2.0;
   RMax_ICF070_1 = 70.0/2.0;
   L_ICF070_1    = 13.0;

   // Frange ICF070(2)
   RMin_ICF070_2 = 20.0/2.0;
   RMax_ICF070_2 = 70.0/2.0;
   L_ICF070_2    = 13.0;

   // Frange V85
   RMin_V85 = 34.0/2.0;
   RMax_V85 = 130.0/2.0;
   L_V85    = 13.0;

   Total_Beam_Length_H  = 0.0;
   Total_Beam_Length_H2 = 0.0;
   Total_Beam_Length_V  = 0.0;

   Total_Beam_Length_H = position_of_beam_start[0];
   // No.1 EGUN  G4Tube
   L_EGUN = 209.0;

    RMin_EGUN = 197.0/2.0;
    RMax_EGUN = RMax_ICF203;
    L_EGUN    = 231.0 - ( L_ICF203 + L_ICF203 );

    RMin_EGUN_FF = 197.0/2.0;
    RMax_EGUN_FF = RMax_ICF203;
    L_EGUN_FF    = L_ICF203; 

    RMin_EGUN_BF = 197.0/2.0;
    RMax_EGUN_BF = RMax_ICF203;
    L_EGUN_BF    = L_ICF203;

    position_EGUN_FF[0] = Total_Beam_Length_H + L_EGUN_FF/2.0;    
    position_EGUN_FF[1] = position_of_beam_start[1];
    position_EGUN_FF[2] = position_of_beam_start[2];

    position_EGUN[0] = Total_Beam_Length_H + L_EGUN_FF + L_EGUN/2.0;
    position_EGUN[1] = position_EGUN_FF[1];
    position_EGUN[2] = position_EGUN_FF[2]; 

    position_EGUN_BF[0] = Total_Beam_Length_H + L_EGUN_FF + L_EGUN + L_EGUN_BF/2.0;
    position_EGUN_BF[1] = position_EGUN[1];
    position_EGUN_BF[2] = position_EGUN[2];

    rotation_EGUN[0] =  0.0;
    rotation_EGUN[1] = 90.0;
    rotation_EGUN[2] =  0.0;

    // Blank Frange            Fragne V85-Blank
    RMin_EGUNBlank = 0.0;  
    RMax_EGUNBlank = RMax_EGUN_FF;
    L_EGUNBlank    = L_ICF203;
    position_EGUNBlank[0] = Total_Beam_Length_H - L_EGUNBlank/2.0;
    position_EGUNBlank[1] = position_of_beam_start[1];
    position_EGUNBlank[2] = position_of_beam_start[2];

    rotation_EGUNBlank[0] =  0.0;
    rotation_EGUNBlank[1] = 90.0;
    rotation_EGUNBlank[2] =  0.0;

    // Anode 
    RMin_Anode = 12.0/2.0;  
    RMax_Anode = 197.0/2.0;  
    L_Anode    = 10.0;
    position_Anode[0] = Total_Beam_Length_H + L_EGUN_FF + L_EGUN - L_Anode/2.0;
    position_Anode[1] = position_of_beam_start[1];
    position_Anode[2] = position_of_beam_start[2];

    rotation_Anode[0] =  0.0;
    rotation_Anode[1] = 90.0;
    rotation_Anode[2] =  0.0;

    Total_Beam_Length_H += L_EGUN + ( L_ICF203*2.0 );
     
   // No.2 ML         G4Tube[0]     Fragne ICF203-ICF070 
    RMin_ML = 29.0/2.0;
    RMax_ML = 35.0/2.0;
    L_ML    = 231.0 - ( L_ICF203 + L_ICF070_2 );

    RMin_ML_FF = RMin_ICF203;
    RMax_ML_FF = RMax_ICF203;
    L_ML_FF    = L_ICF203; 

    RMin_ML_BF = RMin_ICF070_2;
    RMax_ML_BF = RMax_ICF070_2;
    L_ML_BF    = L_ICF070_2;

    position_ML_FF[0] = Total_Beam_Length_H + L_ICF203/2.0;    
    position_ML_FF[1] = position_of_beam_start[1];
    position_ML_FF[2] = position_of_beam_start[2];

    position_ML[0] = Total_Beam_Length_H + L_ICF203 + L_ML/2.0;
    position_ML[1] = position_ML_FF[1];
    position_ML[2] = position_ML_FF[2]; 

    position_ML_BF[0] = Total_Beam_Length_H + L_ICF203 + L_ML + L_ICF070_2/2.0;
    position_ML_BF[1] = position_ML[1];
    position_ML_BF[2] = position_ML[2];

    rotation_ML[0] =  0.0;
    rotation_ML[1] = 90.0;
    rotation_ML[2] =  0.0;

    Total_Beam_Length_H += L_ML + ( L_ICF203 + L_ICF070_2 );      

   // No.3 ICF070-GV  G4Tube        Fragne ICF070-ICF070 
    RMin_ICF070GV = 38.0/2.0;
    RMax_ICF070GV = 70.0/2.0;
    L_ICF070GV    = 35.0 - ( L_ICF070_1 + L_ICF070_1 );

    RMin_ICF070GV_FF = RMin_ICF070_1; 
    RMax_ICF070GV_FF = RMax_ICF070_1;
    L_ICF070GV_FF    = L_ICF070_1;

    RMin_ICF070GV_BF = RMin_ICF070_1;
    RMax_ICF070GV_BF = RMax_ICF070_1;
    L_ICF070GV_BF    = L_ICF070_1;

    position_ICF070GV_FF[0] = Total_Beam_Length_H + L_ICF070_1/2.0;
    position_ICF070GV_FF[1] = position_of_beam_start[1]; 
    position_ICF070GV_FF[2] = position_of_beam_start[2];

    position_ICF070GV[0] = Total_Beam_Length_H + L_ICF070_1 + L_ICF070GV/2.0;
    position_ICF070GV[1] = position_ICF070GV_FF[1];
    position_ICF070GV[2] = position_ICF070GV_FF[2];

    position_ICF070GV_BF[0] = Total_Beam_Length_H + L_ICF070_1 + L_ICF070GV +  L_ICF070_1/2.0;
    position_ICF070GV_BF[1] = position_ICF070GV[1];
    position_ICF070GV_BF[2] = position_ICF070GV[2];

    rotation_ICF070GV[0] =  0.0;
    rotation_ICF070GV[1] = 90.0;
    rotation_ICF070GV[2] =  0.0;

    Total_Beam_Length_H += L_ICF070GV + ( L_ICF070_1 + L_ICF070_1 );

   // No.4 ICF070-V85 duct G4Tube   Fragne ICF070-V85 
    RMin_ICF070V85Duct = 39.4/2.0;
    RMax_ICF070V85Duct = 42.7/2.0;
    L_ICF070V85Duct    = 80.0 - ( L_ICF070_1 + L_V85 );   

    RMin_ICF070V85Duct_FF = RMin_ICF070_1;
    RMax_ICF070V85Duct_FF = RMax_ICF070_1;
    L_ICF070V85Duct_FF    = L_ICF070_1;

    RMin_ICF070V85Duct_BF = RMin_V85;
    RMax_ICF070V85Duct_BF = RMax_V85;
    L_ICF070V85Duct_BF    = L_V85;

    position_ICF070V85Duct_FF[0] = Total_Beam_Length_H + L_ICF070_1/2.0;
    position_ICF070V85Duct_FF[1] = position_of_beam_start[1];
    position_ICF070V85Duct_FF[2] = position_of_beam_start[2];

    position_ICF070V85Duct[0] = Total_Beam_Length_H + L_ICF070_1 + L_ICF070V85Duct/2.0;
    position_ICF070V85Duct[1] = position_ICF070V85Duct_FF[1];
    position_ICF070V85Duct[2] = position_ICF070V85Duct_FF[2];

    position_ICF070V85Duct_BF[0] = Total_Beam_Length_H + L_ICF070_1 + L_ICF070V85Duct + L_V85/2.0;
    position_ICF070V85Duct_BF[1] = position_ICF070V85Duct[1];
    position_ICF070V85Duct_BF[2] = position_ICF070V85Duct[2];

    rotation_ICF070V85Duct[0] =  0.0;
    rotation_ICF070V85Duct[1] = 90.0;
    rotation_ICF070V85Duct[2] =  0.0;

    Total_Beam_Length_H += L_ICF070V85Duct + ( L_ICF070_1 + L_V85 );

   // No.5 CoreMonitor1 G4Tube      Fragne V85-V85
    RMin_CM1 = 60.5/2.0;
    RMax_CM1 = 66.5/2.0;
    L_CM1    = 160.0 - ( L_V85*2.0 );

    RMin_CM1_FF = RMin_V85;
    RMax_CM1_FF = RMax_V85;
    L_CM1_FF    = L_V85;

    RMin_CM1_BF = RMin_V85;
    RMax_CM1_BF = RMax_V85;
    L_CM1_BF    = L_V85;

    position_CM1_FF[0] = Total_Beam_Length_H + L_V85/2.0;
    position_CM1_FF[1] = position_of_beam_start[1];
    position_CM1_FF[2] = position_of_beam_start[2];

    position_CM1[0] = Total_Beam_Length_H + L_V85 + L_CM1/2.0;     
    position_CM1[1] = position_CM1_FF[1];
    position_CM1[2] = position_CM1_FF[2];

    position_CM1_BF[0] = Total_Beam_Length_H + L_V85 + L_CM1 + L_V85/2.0;
    position_CM1_BF[1] = position_CM1[1];
    position_CM1_BF[2] = position_CM1[2];

    rotation_CM1[0] =  0.0;
    rotation_CM1[1] = 90.0;
    rotation_CM1[2] =  0.0;

    Total_Beam_Length_H += L_CM1 + ( L_V85*2.0 );

   // No.6 PB+B-Tube G4Tube         Fragne V85-V85 
    RMin_PBBTube = 20.0/2.0;
    RMax_PBBTube = 130.0/2.0;
    L_PBBTube    = 1438.1 - ( L_V85*2.0 );

    RMin_PBBTube_FF = RMin_V85;
    RMax_PBBTube_FF = RMax_V85;
    L_PBBTube_FF    = L_V85;

    RMin_PBBTube_BF = RMin_V85;
    RMax_PBBTube_BF = RMax_V85;
    L_PBBTube_BF    = L_V85;

    position_PBBTube_FF[0] = Total_Beam_Length_H + L_V85/2.0;
    position_PBBTube_FF[1] = position_of_beam_start[1];
    position_PBBTube_FF[2] = position_of_beam_start[2];

    position_PBBTube[0] = Total_Beam_Length_H + L_V85 + L_PBBTube/2.0;     
    position_PBBTube[1] = position_PBBTube_FF[1];
    position_PBBTube[2] = position_PBBTube_FF[2];

    position_PBBTube_BF[0] = Total_Beam_Length_H + L_V85 + L_PBBTube + L_V85/2.0;
    position_PBBTube_BF[1] = position_PBBTube[1];
    position_PBBTube_BF[2] = position_PBBTube[2];

    rotation_PBBTube[0] =  0.0;
    rotation_PBBTube[1] = 90.0;
    rotation_PBBTube[2] =  0.0;

    Total_Beam_Length_H += L_PBBTube + ( L_V85*2.0 );

    // No.7 CoreMonitor2 G4Tube      Fragne V85-V85
    RMin_CM2 = 60.5/2.0;
    RMax_CM2 = 66.5/2.0;
    L_CM2    = 198.2 - ( L_V85*2.0 );

    RMin_CM2_FF = RMin_V85;
    RMax_CM2_FF = RMax_V85;
    L_CM2_FF     = L_V85;

    RMin_CM2_BF = RMin_V85;
    RMax_CM2_BF = RMax_V85;
    L_CM2_BF     = L_V85;

    position_CM2_FF[0] = Total_Beam_Length_H + L_V85/2.0;
    position_CM2_FF[1] = position_of_beam_start[1];
    position_CM2_FF[2] = position_of_beam_start[2];

    position_CM2[0] = Total_Beam_Length_H + L_V85 + L_CM2/2.0;     
    position_CM2[1] = position_CM2_FF[1];
    position_CM2[2] = position_CM2_FF[2];

    position_CM2_BF[0] = Total_Beam_Length_H + L_V85 + L_CM2 + L_V85/2.0;
    position_CM2_BF[1] = position_CM2[1];
    position_CM2_BF[2] = position_CM2[2];

    rotation_CM2[0] =  0.0;
    rotation_CM2[1] = 90.0;
    rotation_CM2[2] =  0.0;

    Total_Beam_Length_H += L_CM2 + ( L_V85*2.0 );

    // No.8 2m-Tube G4Tube           Fragne V85-V85 
    RMin_2mTube = 20.0/2.0;
    RMax_2mTube = 130.0/2.0;
    L_2mTube    = 2037.5 - ( L_V85*2.0 );

    RMin_2mTube_FF = RMin_V85;
    RMax_2mTube_FF = RMax_V85;
    L_2mTube_FF    = L_V85;

    RMin_2mTube_BF = RMin_V85;
    RMax_2mTube_BF = RMax_V85;
    L_2mTube_BF    = L_V85;

    position_2mTube_FF[0] = Total_Beam_Length_H + L_V85/2.0;
    position_2mTube_FF[1] = position_of_beam_start[1];
    position_2mTube_FF[2] = position_of_beam_start[2];

    position_2mTube[0] = Total_Beam_Length_H + L_V85 + L_2mTube/2.0;     
    position_2mTube[1] = position_2mTube_FF[1];
    position_2mTube[2] = position_2mTube_FF[2];

    position_2mTube_BF[0] = Total_Beam_Length_H + L_V85 + L_2mTube + L_V85/2.0;
    position_2mTube_BF[1] = position_2mTube[1];
    position_2mTube_BF[2] = position_2mTube[2];

    rotation_2mTube[0] =  0.0;
    rotation_2mTube[1] = 90.0;
    rotation_2mTube[2] =  0.0;

    Total_Beam_Length_H += L_2mTube + ( L_V85*2.0 );

   // No.9 Drift Duct G4Tube        Fragne V85-V85
    RMin_V85V85Duct = 36.5/2.0;
    RMax_V85V85Duct = 42.7/2.0;
    L_V85V85Duct    = 124.0 - ( L_V85*2.0 );

    RMin_V85V85Duct_FF = RMin_V85;
    RMax_V85V85Duct_FF = RMax_V85;
    L_V85V85Duct_FF    = L_V85;

    RMin_V85V85Duct_BF = RMin_V85;
    RMax_V85V85Duct_BF = RMax_V85;
    L_V85V85Duct_BF    = L_V85;

    position_V85V85Duct_FF[0] = Total_Beam_Length_H + L_V85/2.0;
    position_V85V85Duct_FF[1] = position_of_beam_start[1];
    position_V85V85Duct_FF[2] = position_of_beam_start[2];

    position_V85V85Duct[0] = Total_Beam_Length_H + L_V85 + L_V85V85Duct/2.0;     
    position_V85V85Duct[1] = position_V85V85Duct_FF[1];
    position_V85V85Duct[2] = position_V85V85Duct_FF[2];

    position_V85V85Duct_BF[0] = Total_Beam_Length_H + L_V85 + L_V85V85Duct + L_V85/2.0;
    position_V85V85Duct_BF[1] = position_V85V85Duct[1];
    position_V85V85Duct_BF[2] = position_V85V85Duct[2];

    rotation_V85V85Duct[0] =  0.0;
    rotation_V85V85Duct[1] = 90.0;
    rotation_V85V85Duct[2] =  0.0;

    Total_Beam_Length_H += L_V85V85Duct + ( L_V85*2.0 );

   // No.10 V85-GV G4Tube           Fragne V85-V85
    RMin_V85GV = 64.0/2.0;
    RMax_V85GV = 130.0/2.0;
    L_V85GV    = 60.0 - ( L_V85*2.0 );

    RMin_V85GV_FF = RMin_V85;
    RMax_V85GV_FF = RMax_V85; 
    L_V85GV_FF    = L_V85;

    RMin_V85GV_BF = RMin_V85;
    RMax_V85GV_BF = RMax_V85; 
    L_V85GV_BF    = L_V85;

    position_V85GV_FF[0] = Total_Beam_Length_H + L_V85/2.0;
    position_V85GV_FF[1] = position_of_beam_start[1];
    position_V85GV_FF[2] = position_of_beam_start[2];

    position_V85GV[0] = Total_Beam_Length_H + L_V85 + L_V85GV/2.0;     
    position_V85GV[1] = position_V85GV_FF[1];
    position_V85GV[2] = position_V85GV_FF[2];

    position_V85GV_BF[0] = Total_Beam_Length_H + L_V85 + L_V85GV + L_V85/2.0;
    position_V85GV_BF[1] = position_V85GV[1];
    position_V85GV_BF[2] = position_V85GV[2];

    rotation_V85GV[0] =  0.0;
    rotation_V85GV[1] = 90.0;
    rotation_V85GV[2] =  0.0;

    Total_Beam_Length_H += L_V85GV + ( L_V85*2.0 );

   // No.11 QM G4Tube               Fragne V85-V85
    RMin_QM = 29.0/2.0;
    RMax_QM = 35.0/2.0;
    L_QM    = 672.5 - ( L_V85*2.0 );

    RMin_QM_FF = RMin_V85;
    RMax_QM_FF = RMax_V85;
    L_QM_FF    = L_V85;

    RMin_QM_BF = RMin_V85;
    RMax_QM_BF = RMax_V85;
    L_QM_BF    = L_V85;

    position_QM_FF[0] = Total_Beam_Length_H + L_V85/2.0;
    position_QM_FF[1] = position_of_beam_start[1];
    position_QM_FF[2] = position_of_beam_start[2];

    position_QM[0] = Total_Beam_Length_H + L_V85 + L_QM/2.0;
    position_QM[1] = position_QM_FF[1];
    position_QM[2] = position_QM_FF[2];

    position_QM_BF[0] = Total_Beam_Length_H + L_V85 + L_QM + L_V85/2.0;
    position_QM_BF[1] = position_QM[1];
    position_QM_BF[2] = position_QM[2];

    rotation_QM[0] =  0.0;
    rotation_QM[1] = 90.0;
    rotation_QM[2] =  0.0;

    Total_Beam_Length_H += L_QM + ( L_V85*2.0 );

   // No.12 ScreenMonitor1  Fragne V85-V85
    L_SM1 = 160.0;

    // parts-0 ( Frange )
    RMin_SM1_FF = RMin_V85;
    RMax_SM1_FF = RMax_V85;
    L_SM1_FF    = L_V85;

    RMin_SM1_BF = RMin_V85;
    RMax_SM1_BF = RMax_V85;
    L_SM1_BF    = L_V85;

    position_SM1_FF[0] = Total_Beam_Length_H + L_V85/2.0;
    position_SM1_FF[1] = position_of_beam_start[1];
    position_SM1_FF[2] = position_of_beam_start[2];

    position_SM1_BF[0] = Total_Beam_Length_H + L_SM1 - L_V85/2.0;
    position_SM1_BF[1] = position_of_beam_start[1];
    position_SM1_BF[2] = position_of_beam_start[2];

    rotation_SM1_F[0] =  0.0;
    rotation_SM1_F[1] = 90.0;
    rotation_SM1_F[2] =  0.0;
            
    // parts-1
    Rmin_SM1_parts1 = 28.0/2.0;
    Rmax_SM1_parts1 = 34.0/2.0;
    L_SM1_parts1    = L_SM1 - ( L_V85*2.0 );
    
    position_SM1_parts1[0] = Total_Beam_Length_H + L_V85 + L_SM1_parts1/2.0;
    position_SM1_parts1[1] = position_of_beam_start[1]; 
    position_SM1_parts1[2] = position_of_beam_start[2];

    // parts-2
    Rmin_SM1_parts2 = 0.0;
    Rmax_SM1_parts2 = 54.0/2.0;
    L_SM1_parts2    = 100.0;
    
    position_SM1_parts2[0] = position_SM1_parts1[0];
    position_SM1_parts2[1] = position_SM1_parts1[1];
    position_SM1_parts2[2] = position_SM1_parts1[2];
   
    // parts-3
    Rmin_SM1_parts3 = 54.0/2.0;
    Rmax_SM1_parts3 = 60.0/2.0;
    L_SM1_parts3    = 100.0;

    position_SM1_parts3[0] = position_SM1_parts1[0];
    position_SM1_parts3[1] = position_SM1_parts1[1];
    position_SM1_parts3[2] = position_SM1_parts1[2];

    rotation_SM1_comp[0] = 0.0;
    rotation_SM1_comp[1] = 90.0;
    rotation_SM1_comp[2] = 0.0;

    rotation_SM1_duct[0] =  0.0;
    rotation_SM1_duct[1] = -90.0;
    rotation_SM1_duct[2] =  0.0;

    rotation_SM1_body[0] = 0.0;
    rotation_SM1_body[1] = 0.0;
    rotation_SM1_body[2] = 0.0;

    Total_Beam_Length_H += L_SM1; 

      // No.12.5 Alumina AF995R   
  
   // No.13 CoreMonitor3 G4Tube     Fragne V85-V85
    RMin_CM3 = 60.5/2.0;
    RMax_CM3 = 66.5/2.0;
    L_CM3    = 200.0 - ( L_V85*2.0 );
    RMin_CM3_FF = RMin_V85;
    RMax_CM3_FF = RMax_V85;
    L_CM3_FF    = L_V85;

    RMin_CM3_BF = RMin_V85;
    RMax_CM3_BF = RMax_V85;
    L_CM3_BF    = L_V85;

    position_CM3_FF[0] = Total_Beam_Length_H + L_V85/2.0;
    position_CM3_FF[1] = position_of_beam_start[1];
    position_CM3_FF[2] = position_of_beam_start[2];

    position_CM3[0] = Total_Beam_Length_H + L_V85 + L_CM3/2.0;
    position_CM3[1] = position_CM3_FF[1];
    position_CM3[2] = position_CM3_FF[2];

    position_CM3_BF[0] = Total_Beam_Length_H + L_V85 + L_CM3 + L_V85/2.0;
    position_CM3_BF[1] = position_CM3[1];
    position_CM3_BF[2] = position_CM3[2];

    rotation_CM3[0] =  0.0;
    rotation_CM3[1] = 90.0;
    rotation_CM3[2] =  0.0;

    Total_Beam_Length_H += L_CM3 + ( L_V85*2.0 );
  
   // No.14 BM G4                   Fragne V85-V85-V85 
    L1_BM = 377.0;
    LH_BM = 604.0;
    LV_BM = 377.0;
 
    // No.14.5 Cu-Collimeter G4
    L_CuCollimeter = 22.0;
    W_CuCollimeter = 16.0;
    H_CuCollimeter = 64.0;
    
    RMin_CuCollimeter_hole = 0.0;
    RMax_CuCollimeter_hole = 10.0/2.0;   // <-- Real value but it can be changed by user.
    L_CuCollimeter_hole = L_CuCollimeter;

    size_of_CuCollimeter[0] = L_CuCollimeter;
    size_of_CuCollimeter[1] = W_CuCollimeter;
    size_of_CuCollimeter[2] = H_CuCollimeter;
    
    position_CuCollimeter[0] = Total_Beam_Length_H + L_CuCollimeter/2.0;
    position_CuCollimeter[1] = position_of_beam_start[1];
    position_CuCollimeter[2] = position_of_beam_start[2];

    cout << "Cu-Collimeter Poistion = " 
         <<  position_CuCollimeter[0]  << " " 
	 <<  position_CuCollimeter[1]  << " "
         <<  position_CuCollimeter[2]  << endl;   

    if( fabs(alignment_position_CuCollimeter_hole[2])>H_CuCollimeter/2.0 )
      { 
        alignment_position_CuCollimeter_hole[2] = 0.0; 
        fprintf(stderr, "The shift of CuCollimeter hole position is over the maximum value ...\n");
        fprintf(stderr, "The shift becomes to be 0...\n");
       }
 
    for( int i=0; i<3; i++ ){
         position_CuCollimeter_hole[i] = position_CuCollimeter[i] + alignment_position_CuCollimeter_hole[i];
    }
    rotation_CuCollimeter[0] =  0.0;
    rotation_CuCollimeter[1] =  0.0;
    rotation_CuCollimeter[2] =  0.0;

    rotation_CuCollimeter_hole[0] =   0.0;
    rotation_CuCollimeter_hole[1] =  90.0;
    rotation_CuCollimeter_hole[2] =   0.0;

    // part-0 BM-Duct Frange 
    // Three Franges
    // 0 : Insident frange
    // 1 : Straight out frange
    // 2 : Bending out frange
    Rmin_V85_BM=0;
    Rmax_V85_BM=130.0/2.0;
    L_V85_BM=12.0;
    
    size_of_V85_BM_IN[0]=64.0;
    size_of_V85_BM_IN[1]=16.0;
    size_of_V85_BM_IN[2]=L_V85_BM;

    position_BM_F[0][0]=Total_Beam_Length_H+L_V85_BM/2.0;
    position_BM_F[0][1]=position_of_beam_start[1];
    position_BM_F[0][2]=position_of_beam_start[2];

    rotation_BM_F[0][0]=0.0;
    rotation_BM_F[0][1]=90.0;
    rotation_BM_F[0][2]=0.0;

    position_BM_F[1][0]=Total_Beam_Length_H + LH_BM - L_V85_BM/2.0;
    position_BM_F[1][1]=position_of_beam_start[1];
    position_BM_F[1][2]=position_of_beam_start[2];

    rotation_BM_F[1][0]=0.0;
    rotation_BM_F[1][1]=90.0;
    rotation_BM_F[1][2]=0.0;

    position_BM_F[2][0]=Total_Beam_Length_H + L1_BM;
    position_BM_F[2][1]=position_of_beam_start[1];
    position_BM_F[2][2]=position_of_beam_start[2] + LV_BM - L_V85_BM/2.0;

    rotation_BM_F[2][0]=0.0;
    rotation_BM_F[2][1]=0.0;
    rotation_BM_F[2][2]=0.0;

    // parts-1 : Straight Duct ( Insident part )
    size_of_BM_parts_1_outer[0] = 165.0;
    size_of_BM_parts_1_outer[1] = 22.0;
    size_of_BM_parts_1_outer[2] = 70.0;

    size_of_BM_parts_1_inner[0] = size_of_BM_parts_1_outer[0];
    size_of_BM_parts_1_inner[1] = 16.0;
    size_of_BM_parts_1_inner[2] = 64.0;

    position_of_BM_parts_1_outer[0] = Total_Beam_Length_H + L_V85_BM 
                                      + size_of_BM_parts_1_outer[0]/2.0;
    position_of_BM_parts_1_outer[1] = position_of_beam_start[1];
    position_of_BM_parts_1_outer[2] = position_of_beam_start[2];
   
    for( int i=0; i<3; i++ ) position_of_BM_parts_1_inner[i]=position_of_BM_parts_1_outer[i]; 

    // parts-1 : Straight Duct vacuum (  Straight entrance part )
    size_of_BM_parts_1_vac[0] = 155.0;  // 165 + 12 - 22 = 155.0
    size_of_BM_parts_1_vac[1] = size_of_BM_parts_1_inner[1];
    size_of_BM_parts_1_vac[2] = size_of_BM_parts_1_inner[2];

    position_of_BM_parts_1_vac[0] = Total_Beam_Length_H + L_CuCollimeter 
                                     + size_of_BM_parts_1_vac[0]/2.0;
    position_of_BM_parts_1_vac[1] = position_of_BM_parts_1_outer[1];   
    position_of_BM_parts_1_vac[2] = position_of_BM_parts_1_outer[2];

    // parts-2 : Straight Duct ( Straighte exit part )
    size_of_BM_parts_2_outer[0] = 130.0;
    size_of_BM_parts_2_outer[1] = 22.0;
    size_of_BM_parts_2_outer[2] = 70.0;
    
    size_of_BM_parts_2_inner[0] = size_of_BM_parts_2_outer[0];
    size_of_BM_parts_2_inner[1] = 16.0;
    size_of_BM_parts_2_inner[2] = 64.0;

    position_of_BM_parts_2_outer[0] = Total_Beam_Length_H + LH_BM 
                                      - L_V85_BM
                                      - size_of_BM_parts_2_outer[0]/2.0;
    position_of_BM_parts_2_outer[1] = position_of_beam_start[1];
    position_of_BM_parts_2_outer[2] = position_of_beam_start[2];

    for( int i=0; i<3; i++ ) position_of_BM_parts_2_inner[i]=position_of_BM_parts_2_outer[i];

    // parts-2 : Straight Duct vacuum
    for( int i=0; i<3; i++ ){
      size_of_BM_parts_2_vac[i] = size_of_BM_parts_2_inner[i];
      position_of_BM_parts_2_vac[i] = position_of_BM_parts_2_inner[i];
    } 
 
    // parts-3 : Straight Duct ( Bending exit part )
    size_of_BM_parts_3_outer[0] = 70.0;
    size_of_BM_parts_3_outer[1] = 22.0;
    size_of_BM_parts_3_outer[2] = 165.0;

    size_of_BM_parts_3_inner[0] = 64.0;
    size_of_BM_parts_3_inner[1] = 16.0;
    size_of_BM_parts_3_inner[2] = size_of_BM_parts_3_outer[2];

    position_of_BM_parts_3_outer[0] = Total_Beam_Length_H + L1_BM;
    position_of_BM_parts_3_outer[1] = position_of_beam_start[1];
    position_of_BM_parts_3_outer[2] = position_of_beam_start[2] + LV_BM
                                       - L_V85_BM
                                       - size_of_BM_parts_3_outer[2]/2.0;
    for( int i=0; i<3; i++ ) position_of_BM_parts_3_inner[i] = position_of_BM_parts_3_outer[i];

    // parts-3 : Straight Duct vacuum
    for( int i=0; i<3; i++ ){
      size_of_BM_parts_3_vac[i] = size_of_BM_parts_3_inner[i];
      position_of_BM_parts_3_vac[i] = position_of_BM_parts_3_inner[i];
    }
 
    // parts-4(1) : Main Body ( Side Board left )
    size_of_BM_parts_4_1[0] = 285.0;
    size_of_BM_parts_4_1[1] = 3.0;
    size_of_BM_parts_4_1[2] = 235.0; 
 
    position_of_BM_parts_4_1[0] = Total_Beam_Length_H 
                                   + L_V85_BM 
                                   + size_of_BM_parts_1_outer[0]
                                   + size_of_BM_parts_4_1[0]/2.0;
    position_of_BM_parts_4_1[1] = position_of_beam_start[1] 
                                  - 8.0 
                                  - size_of_BM_parts_4_1[1]/2.0;
    position_of_BM_parts_4_1[2] = position_of_beam_start[2] 
                                  + 82.5;  // 82.5 = 200.0 - 235.0/2.0 = 235.0/2. 0- 35  
   
    // parts-4(2) : Main Body ( Side Board right )
    size_of_BM_parts_4_2[0] = size_of_BM_parts_4_1[0];
    size_of_BM_parts_4_2[1] = size_of_BM_parts_4_1[1];
    size_of_BM_parts_4_2[2] = size_of_BM_parts_4_1[2];

    position_of_BM_parts_4_2[0] = position_of_BM_parts_4_1[0];
    position_of_BM_parts_4_2[1] = position_of_beam_start[1]
                                  + 8.0 
                                  + size_of_BM_parts_4_1[1]/2.0; 
    position_of_BM_parts_4_2[2] = position_of_BM_parts_4_1[2];

    // parts-4(3) : Main Body ( Side Board Cut Volume 1 )
    Rmin_BM_parts_4_3 = 0.0;
    Rmax_BM_parts_4_3 = 165.0;
    L_BM_parts_4_3    = 30.0;
    position_of_BM_parts_4_3[0] =  Total_Beam_Length_H
                                   + L_V85_BM
                                   + size_of_BM_parts_1_outer[0];
    position_of_BM_parts_4_3[1] = position_of_beam_start[1];
    position_of_BM_parts_4_3[2] = position_of_beam_start[2]
                                 + LV_BM
                                 - L_V85_BM
                                 - size_of_BM_parts_3_outer[2];
    rotation_of_BM_parts_4_3[0] = 90.0;
    rotation_of_BM_parts_4_3[1] = 0.0;
    rotation_of_BM_parts_4_3[2] = 0.0;    

    // parts-4(4) : Main Body ( Side Board Cut Volume 2 )
    size_of_BM_parts_4_4[0] = 50.0+10.0;
    size_of_BM_parts_4_4[1] = 22.0+18.0;
    size_of_BM_parts_4_4[2] = 115.0+5.0;
    position_of_BM_parts_4_4[0] =  Total_Beam_Length_H
                                   + L_V85_BM
                                   + size_of_BM_parts_1_outer[0]
                                   + size_of_BM_parts_4_1[0]
                                   - size_of_BM_parts_4_4[0]/2.0
                                   + ( size_of_BM_parts_4_4[0]-50.0 );
    position_of_BM_parts_4_4[1] = position_of_beam_start[1];
    position_of_BM_parts_4_4[2] = position_of_BM_parts_4_3[2]
                                  - size_of_BM_parts_4_4[2]/2.0
                                  + ( size_of_BM_parts_4_4[2]-115.0 );

    // parts-4(5) : Main Body ( Side Board Cut Volume 3 )
    Rmin_BM_parts_4_5 = 0.0;
    Rmax_BM_parts_4_5 = 50.0;
    L_BM_parts_4_5 = 30.0;
    position_of_BM_parts_4_5[0] = Total_Beam_Length_H 
                                  + L_V85_BM
                                  + size_of_BM_parts_1_outer[0]
                                  + size_of_BM_parts_4_1[0];
    position_of_BM_parts_4_5[1] = position_of_beam_start[1];
    position_of_BM_parts_4_5[2] = position_of_beam_start[2]
                                  + LV_BM
                                  - L_V85_BM
                                  - size_of_BM_parts_3_outer[2]
                                  - size_of_BM_parts_4_4[2];

    rotation_of_BM_parts_4_5[0] = 90.0;
    rotation_of_BM_parts_4_5[1] = 0.0;
    rotation_of_BM_parts_4_5[2] = 0.0;    

    // parts-4(6) : Main Body ( Side Board  )
    Rmin_BM_parts_4_6 = 165.0;
    Rmax_BM_parts_4_6 = 168.0;
    L_BM_parts_4_6    = 16.0;
    Phi_start_BM_parts_4_6 = 270.0;
    Phi_end_BM_parts_4_6  = 90.0;
    for( int i=0; i<3; i++ ) position_of_BM_parts_4_6[i] = position_of_BM_parts_4_3[i];
    for( int i=0; i<3; i++ ) rotation_of_BM_parts_4_6[i] = rotation_of_BM_parts_4_3[i];

    rotation_of_BM_parts_4_6[0] = 90.0;
    rotation_of_BM_parts_4_6[1] = 0.0;
    rotation_of_BM_parts_4_6[2] = 0.0;    
    
    // parts-4(7) : Main Body ( Side Board  )
    size_of_BM_parts_4_7[0] = 3.0;
    size_of_BM_parts_4_7[1] = 16.0;
    size_of_BM_parts_4_7[2] = 115.0; 
    position_of_BM_parts_4_7[0] = position_of_BM_parts_4_4[0]
                                  - size_of_BM_parts_4_4[0]/2.0
                                  - size_of_BM_parts_4_7[0]/2.0;
    position_of_BM_parts_4_7[1] = position_of_BM_parts_4_4[1];
    position_of_BM_parts_4_7[2] = position_of_BM_parts_4_4[2];

    // parts-4(8) : Main Body ( Side Board  )
    Rmin_BM_parts_4_8 = 50.0;
    Rmax_BM_parts_4_8 = 53.0;
    L_BM_parts_4_8    = 16.0;
    Phi_start_BM_parts_4_8 = 180.0;
    Phi_end_BM_parts_4_8   = 90.0;
    for( int i=0; i<3; i++ ) position_of_BM_parts_4_8[i] = position_of_BM_parts_4_5[i];
    for( int i=0; i<3; i++ ) rotation_of_BM_parts_4_8[i] = rotation_of_BM_parts_4_5[i];
    rotation_of_BM_parts_4_8[0] = 90.0;
    rotation_of_BM_parts_4_8[1] = 0.0;
    rotation_of_BM_parts_4_8[2] = 0.0;    

    // parts-4(9) : Main Body ( Bottom Board  )
    size_of_BM_parts_4_9[0] = 285.0;
    size_of_BM_parts_4_9[1] = 16.0;
    size_of_BM_parts_4_9[2] = 3.0;
    position_of_BM_parts_4_9[0] = Total_Beam_Length_H
                                  + L_V85_BM
                                  + size_of_BM_parts_1_outer[0]
                                  + size_of_BM_parts_4_9[0]/2.0;
    position_of_BM_parts_4_9[1] = position_of_beam_start[1];
    position_of_BM_parts_4_9[2] = position_of_beam_start[2]
                                  - 35.0
                                  + size_of_BM_parts_4_9[2]/2.0;

   //Vacuum region of BM-Duct

    // parts-5(1) : Main Body
    size_of_BM_parts_5_1[0] = 285.0;
    size_of_BM_parts_5_1[1] = 16.0;
    size_of_BM_parts_5_1[2] = 232.0; // 232.0 = 235.0 - 3.0   
    position_of_BM_parts_5_1[0] = position_of_BM_parts_4_1[0];
    position_of_BM_parts_5_1[1] = position_of_beam_start[1];
    position_of_BM_parts_5_1[2] = position_of_beam_start[2]
                                  + 84.0;  // 84.0 = 200.0 - 232.0/2.0 = 232.0/2.0 - 32  
   
    // parts-5(2) : Main Body ( Cut Volume 1 )
    Rmin_BM_parts_5_2 = 0.0;
    Rmax_BM_parts_5_2 = 168.0;
    L_BM_parts_5_2    = 16.0 + 10.0;
    for( int i=0; i<3; i++ ) position_of_BM_parts_5_2[i] = position_of_BM_parts_4_3[i];
    for( int i=0; i<3; i++ ) rotation_of_BM_parts_5_2[i] = rotation_of_BM_parts_4_3[i]; 
    
    // parts-5(3) : Main Body ( Cut Volume 2 )
    size_of_BM_parts_5_3[0] = 53.0 + 10.0;
    size_of_BM_parts_5_3[1] = 16.0 + 10.0;
    size_of_BM_parts_5_3[2] = 115.0 + 10.0;
    position_of_BM_parts_5_3[0] = position_of_BM_parts_5_1[0]
                                  + size_of_BM_parts_5_1[0]/2.0
                                  - size_of_BM_parts_5_3[0]/2.0
                                  + ( size_of_BM_parts_5_3[0]-53.0 ); 
    position_of_BM_parts_5_3[1] = position_of_beam_start[1];
    position_of_BM_parts_5_3[2] = position_of_BM_parts_5_2[2]
                                  - size_of_BM_parts_5_3[2]/2.0
                                  + ( size_of_BM_parts_5_3[2]-115.0 );

    // parts-5(4) : Main Body ( Cut Volume 3 )
    Rmin_BM_parts_5_4 = 0.0;
    Rmax_BM_parts_5_4 = 53.0;
    L_BM_parts_5_4    = 16.0 + 10.0;
    for( int i=0; i<3; i++ ) position_of_BM_parts_5_4[i] = position_of_BM_parts_4_5[i];
    for( int i=0; i<3; i++ ) rotation_of_BM_parts_5_4[i] = rotation_of_BM_parts_4_5[i];

    Total_Beam_Length_H2 = Total_Beam_Length_H + L1_BM; 
    Total_Beam_Length_H += LH_BM; 
    Total_Beam_Length_V = position_of_beam_start[2]+LV_BM; 


      // BM York & Coil
      //  BM York parts-1 & 2 
      size_of_BMYork[0][0] = size_of_BMYork[1][0] = 250.0;
      size_of_BMYork[0][1] = size_of_BMYork[1][1] = 270.0;
      size_of_BMYork[0][2] = size_of_BMYork[1][2] = 250.0;
       
      position_of_BMYork[0][0] = position_of_BMYork[1][0] = Total_Beam_Length_H2 - 75.0;

      position_of_BMYork[0][1] = position_of_beam_start[1] + 12 + size_of_BMYork[0][1]/2.0;
      position_of_BMYork[1][1] = position_of_beam_start[1] - 12 - size_of_BMYork[0][1]/2.0;

      position_of_BMYork[0][2] = position_of_BMYork[1][2] = position_of_beam_start[2] + 75.0;

      //  BM York parts-3 & 4
      size_of_BMYork[2][0] = size_of_BMYork[3][0] = 250.0;
      size_of_BMYork[2][1] = size_of_BMYork[3][1] = 150.0;
      size_of_BMYork[2][2] = size_of_BMYork[3][2] = 200.0;

      position_of_BMYork[2][0] = position_of_BMYork[3][0] = Total_Beam_Length_H2 - 75.0;

      position_of_BMYork[2][1] = position_of_beam_start[1] 
                                 - 12 - 120.0 - size_of_BMYork[2][1]/2.0;
      position_of_BMYork[3][1] = position_of_beam_start[1] 
                                 + 12 + 120.0 + size_of_BMYork[3][1]/2.0;
      position_of_BMYork[2][2] = position_of_BMYork[3][2] 
                               = position_of_BMYork[0][2] - size_of_BMYork[0][2]/2.0
                       	         - size_of_BMYork[2][2]/2.0;
      
      //  BM York parts-5
      size_of_BMYork[4][0] = 250.0;
      size_of_BMYork[4][1] = 564.0;
      size_of_BMYork[4][2] = 75.0;
      position_of_BMYork[4][0] = Total_Beam_Length_H2 -75.0;
      position_of_BMYork[4][1] = position_of_beam_start[1];
      position_of_BMYork[4][2] = position_of_BMYork[2][2] - size_of_BMYork[2][2]/2.0
                        	 - size_of_BMYork[4][2]/2.0; 

      // Coil
      size_of_BMCoil_Outer[0] = 510.0;
      size_of_BMCoil_Outer[1] = 110.0;
      size_of_BMCoil_Outer[2] = 495.0;

      size_of_BMCoil_Inner[0] = 510.0 - 250.0;
      size_of_BMCoil_Inner[1] = 110.0 + 40.0;
      size_of_BMCoil_Inner[2] = 495.0 - 250.0;

      position_of_BMCoil[0][0] = position_of_BMCoil[1][0] = Total_Beam_Length_H2 -75.0;
      position_of_BMCoil[0][1] = position_of_beam_start[1]
	                         - 12.0 - 7.5 - size_of_BMCoil_Outer[1]/2.0;
      position_of_BMCoil[1][1] = position_of_beam_start[1]
                         	 + 12.0 + 7.5 + size_of_BMCoil_Outer[1]/2.0;
      position_of_BMCoil[0][2] = position_of_BMCoil[1][2] = position_of_beam_start[2] + 75.0;

      // alignment of BM York & Coil
      if( alignment_position_BMYork_Coil[0] < -35.0 ){
	alignment_position_BMYork_Coil[0] = 0.;
      }
      if( alignment_position_BMYork_Coil[2] > 35.0 ){
        alignment_position_BMYork_Coil[2] = 0.;
      }
      for( int j=0; j<5; j++ ){  
        for( int i=0; i<3; i++ ){
         position_of_BMYork[j][i]+=alignment_position_BMYork_Coil[i];
	}
      }
        for( int i=0; i<3; i++ ){
         position_of_BMCoil[0][i]+=alignment_position_BMYork_Coil[i];        
       }
 
   // No.15 V85-GV2 G4Tube           Fragne V85-V85
    RMin_V85GV2 = 64.0/2.0;
    RMax_V85GV2 = 130.0/2.0;
    L_V85GV2    = 60.0 - ( L_V85*2.0 );

    RMin_V85GV2_FF = RMin_V85;
    RMax_V85GV2_FF = RMax_V85;
    L_V85GV2_FF    = L_V85;

    RMin_V85GV2_BF = RMin_V85;
    RMax_V85GV2_BF = RMax_V85;
    L_V85GV2_BF    = L_V85;

    position_V85GV2_FF[0] = Total_Beam_Length_H + L_V85/2.0;
    position_V85GV2_FF[1] = position_of_beam_start[1];
    position_V85GV2_FF[2] = position_of_beam_start[2];

    position_V85GV2[0] = Total_Beam_Length_H + L_V85 + L_V85GV2/2.0;
    position_V85GV2[1] = position_V85GV2_FF[1];
    position_V85GV2[2] = position_V85GV2_FF[2];

    position_V85GV2_BF[0] = Total_Beam_Length_H + L_V85 + L_V85GV2 + L_V85/2.0;
    position_V85GV2_BF[1] = position_V85GV2[1];
    position_V85GV2_BF[2] = position_V85GV2[2];

    rotation_V85GV2[0] =  0.0;
    rotation_V85GV2[1] = 90.0;
    rotation_V85GV2[2] =  0.0;

    Total_Beam_Length_H+=L_V85GV2 + ( L_V85*2.0 );

   // No.16 T-Duct G4Tube           Fragne V85-V85
    RMin_TDuct = 64.0/2.0;  // <--- Temporary Value
    RMax_TDuct = 70.0/2.0;  // <--- Temporary Value
    L_TDuct    = 135.5 - ( L_V85*2.0 );
 
    RMin_TDuct_FF = RMin_V85;
    RMax_TDuct_FF = RMax_V85; 
    L_TDuct_FF    = L_V85;

    RMin_TDuct_BF = RMin_V85;
    RMax_TDuct_BF = RMax_V85;
    L_TDuct_BF    = L_V85;

    position_TDuct_FF[0] = Total_Beam_Length_H + L_V85/2.0;
    position_TDuct_FF[1] = position_of_beam_start[1];
    position_TDuct_FF[2] = position_of_beam_start[2];

    position_TDuct[0] = Total_Beam_Length_H + L_V85 + L_TDuct/2.0;
    position_TDuct[1] = position_TDuct_FF[1];
    position_TDuct[2] = position_TDuct_FF[2];

    position_TDuct_BF[0] = Total_Beam_Length_H + L_V85 + L_TDuct + L_V85/2.0;
    position_TDuct_BF[1] = position_TDuct[1];
    position_TDuct_BF[2] = position_TDuct[2];

    rotation_TDuct[0] =  0.0;
    rotation_TDuct[1] = 90.0;
    rotation_TDuct[2] =  0.0;

    Total_Beam_Length_H += L_TDuct + ( L_V85*2.0 );

   // No.17 Blank Frange            Fragne V85-Blank
    RMin_Blank = 0.0;  
    RMax_Blank = 130.0/2.0;  
    L_Blank    = 12.0;
    position_Blank[0] = Total_Beam_Length_H + L_Blank/2.0;
    position_Blank[1] = position_of_beam_start[1];
    position_Blank[2] = position_of_beam_start[2];

    rotation_Blank[0] =  0.0;
    rotation_Blank[1] = 90.0;
    rotation_Blank[2] =  0.0;

    Total_Beam_Length_H += L_Blank;

   //==========================
   // Vertical Beam Line
   //==========================

    /* alignment of vertical beam line                                                                                         
     [0] ... shift along x-axis                                                                                              
     [1] ... shift along y-axis */
    Total_Beam_Length_H2+=alignment_vertical_beamline_shift[0]; 
    position_of_beam_start[1]+=alignment_vertical_beamline_shift[1];

   // No.18 DriftTube G4Tube        Fragne V85-V85 
    RMin_V85V85Duct2 = 60.4/2.0;
    RMax_V85V85Duct2 = 66.4/2.0;
    L_V85V85Duct2    = 100.0 - ( L_V85*2.0 );

    RMin_V85V85Duct2_FF = RMin_V85;
    RMax_V85V85Duct2_FF = RMax_V85;
    L_V85V85Duct2_FF    = L_V85;

    RMin_V85V85Duct2_BF = RMin_V85;
    RMax_V85V85Duct2_BF = RMax_V85;
    L_V85V85Duct2_BF    = L_V85;

    position_V85V85Duct2_FF[0] = Total_Beam_Length_H2;
    position_V85V85Duct2_FF[1] = position_of_beam_start[1];
    position_V85V85Duct2_FF[2] = Total_Beam_Length_V + L_V85/2.0;

    position_V85V85Duct2[0] = Total_Beam_Length_H2;
    position_V85V85Duct2[1] = position_V85V85Duct2_FF[1]; 
    position_V85V85Duct2[2] = Total_Beam_Length_V + L_V85+ L_V85V85Duct2/2.0;

    position_V85V85Duct2_BF[0] = Total_Beam_Length_H2;
    position_V85V85Duct2_BF[1] = position_V85V85Duct2[1];
    position_V85V85Duct2_BF[2] = Total_Beam_Length_V + L_V85 + L_V85V85Duct2 + L_V85/2.0;

    rotation_V85V85Duct2[0] =  0.0;
    rotation_V85V85Duct2[1] =  0.0;
    rotation_V85V85Duct2[2] =  0.0;

    Total_Beam_Length_V += L_V85V85Duct2 + (L_V85*2.0);

   // No.19 Slit                    Fragne V85-V85
    L_SLIT = 148.0;

    // parts-0 ( Frange )
    RMin_SLIT_FF = RMin_V85;
    RMax_SLIT_FF = RMax_V85;
    L_SLIT_FF    = L_V85;

    RMin_SLIT_BF = RMin_V85;
    RMax_SLIT_BF = RMax_V85;
    L_SLIT_BF    = L_V85;

    position_SLIT_FF[0] = Total_Beam_Length_H2;
    position_SLIT_FF[1] = position_of_beam_start[1];
    position_SLIT_FF[2] = Total_Beam_Length_V + L_V85/2.0;

    position_SLIT_BF[0] = Total_Beam_Length_H2;
    position_SLIT_BF[1] = position_of_beam_start[1];
    position_SLIT_BF[2] = Total_Beam_Length_V + L_SLIT - L_V85/2.0; 

    rotation_SLIT_F[0] =  0.0;
    rotation_SLIT_F[1] =  0.0;
    rotation_SLIT_F[2] =  0.0;

    // parts-1 
    Rmin_SLIT_parts1 = 36.5/2.0;
    Rmax_SLIT_parts1 = 42.7/2.0;
    L_SLIT_parts1    = L_SLIT - ( L_V85*2.0 );
    
    position_SLIT_parts1[0] = Total_Beam_Length_H2;
    position_SLIT_parts1[1] = position_of_beam_start[1]; 
    position_SLIT_parts1[2] = Total_Beam_Length_V + L_V85 + L_SLIT_parts1/2.0;

    // parts-2
    Rmin_SLIT_parts2 = 0.0;
    Rmax_SLIT_parts2 = 83.3/2.0;
    L_SLIT_parts2    = 250.0;
    
    position_SLIT_parts2[0] = position_SLIT_parts1[0];
    position_SLIT_parts2[1] = position_SLIT_parts1[1];
    position_SLIT_parts2[2] = position_SLIT_parts1[2];

    // parts-3
    Rmin_SLIT_parts3 = 83.3/2.0;
    Rmax_SLIT_parts3 = 89.5/2.0;
    L_SLIT_parts3    = 250.0;

    position_SLIT_parts3[0] = position_SLIT_parts1[0];
    position_SLIT_parts3[1] = position_SLIT_parts1[1];
    position_SLIT_parts3[2] = position_SLIT_parts1[2];

    rotation_SLIT_comp[0] =  0.0;
    rotation_SLIT_comp[1] = 90.0;
    rotation_SLIT_comp[2] =  0.0;

    rotation_SLIT_duct[0] =  0.0;
    rotation_SLIT_duct[1] =  0.0;
    rotation_SLIT_duct[2] =  0.0;

    rotation_SLIT_body[0] =  0.0;
    rotation_SLIT_body[1] = 90.0;
    rotation_SLIT_body[2] =  0.0;

    // Collimeter
    size_of_collimeter[0] = 50.0;
    size_of_collimeter[1] = 50.0;
    size_of_collimeter[2] = 50.0;
 
    size_of_collimeter_vac[0] = size_of_collimeter[2];
    size_of_collimeter_vac[1] = size_of_collimeter[1];
    size_of_collimeter_vac[2] = size_of_collimeter[0];

    if( slit_width/2.0 + alignment_position_Slit_Collimeter[0] < 0 ){
        alignment_position_Slit_Collimeter[0]=0.0;
	fprintf(stderr, "The shift of CSlit collimator is over the maximum value ...\n");
        fprintf(stderr, "The shift becomes to be 0...\n");
    }
    if( slit_width/2.0 - alignment_position_Slit_Collimeter[1] < 0 ){
      alignment_position_Slit_Collimeter[1]=0.0;
      fprintf(stderr, "The shift of CSlit collimator is over the maximum value ...\n");
      fprintf(stderr, "The shift becomes to be 0...\n");
    }

    position_of_collimeter[0][0] = Total_Beam_Length_H2 + size_of_collimeter[0]/2.0
                                      + alignment_position_Slit_Collimeter[0];    
    position_of_collimeter[0][1] = position_of_beam_start[1];
    position_of_collimeter[0][2] = Total_Beam_Length_V + L_V85 + L_SLIT_parts1/2.0 + 1.0;
    /* The collimators are not installed at the center of slit duct. 
         The position along beam line is center of slit duct + 1.0mm */

    position_of_collimeter[1][0] = Total_Beam_Length_H2 - size_of_collimeter[0]/2.0
                                      + alignment_position_Slit_Collimeter[1];    
    position_of_collimeter[1][1] = position_of_collimeter[0][1];
    position_of_collimeter[1][2] = position_of_collimeter[0][2];

    position_of_collimeter_vac[0][0] = Total_Beam_Length_H2 - 1.0;
    position_of_collimeter_vac[0][1] = position_of_beam_start[1];
    position_of_collimeter_vac[0][2] = position_SLIT_parts2[2] + size_of_collimeter_vac[2]/2.0
                                        + alignment_position_Slit_Collimeter[0];
 
    position_of_collimeter_vac[1][0] = position_of_collimeter_vac[0][0];
    position_of_collimeter_vac[1][1] = position_of_collimeter_vac[0][1];
    position_of_collimeter_vac[1][2] = position_SLIT_parts2[2] - size_of_collimeter_vac[2]/2.0
                                        + alignment_position_Slit_Collimeter[1];

    Total_Beam_Length_V += L_SLIT;

   // No.20 Drift Duct G4Tube       Fragne V85-V85
    RMin_V85V85Duct3 = 60.4/2.0;
    RMax_V85V85Duct3 = 66.4/2.0;
    L_V85V85Duct3    = 200.0 - ( L_V85*2.0 );

    RMin_V85V85Duct3_FF = RMin_V85;
    RMax_V85V85Duct3_FF = RMax_V85;
    L_V85V85Duct3_FF    = L_V85;

    RMin_V85V85Duct3_BF = RMin_V85;
    RMax_V85V85Duct3_BF = RMax_V85;
    L_V85V85Duct3_BF    = L_V85;

    position_V85V85Duct3_FF[0] = Total_Beam_Length_H2;
    position_V85V85Duct3_FF[1] = position_of_beam_start[1];
    position_V85V85Duct3_FF[2] = Total_Beam_Length_V + L_V85/2.0;

    position_V85V85Duct3[0] = Total_Beam_Length_H2;
    position_V85V85Duct3[1] = position_V85V85Duct3_FF[1]; 
    position_V85V85Duct3[2] = Total_Beam_Length_V + L_V85+ L_V85V85Duct3/2.0;

    position_V85V85Duct3_BF[0] = Total_Beam_Length_H2;
    position_V85V85Duct3_BF[1] = position_V85V85Duct3[1];
    position_V85V85Duct3_BF[2] = Total_Beam_Length_V + L_V85 + L_V85V85Duct3 + L_V85/2.0;

    rotation_V85V85Duct3[0] =  0.0;
    rotation_V85V85Duct3[1] =  0.0;
    rotation_V85V85Duct3[2] =  0.0;

    Total_Beam_Length_V += L_V85V85Duct3 + (L_V85*2.0);

   // No.21 ScreenMonitor2          Fragne V85-V85
    L_SM2 = 160.0;

    // parts-0 ( Frange )
    RMin_SM2_FF = RMin_V85;
    RMax_SM2_FF = RMax_V85;
    L_SM2_FF    = L_V85;

    RMin_SM2_BF = RMin_V85;
    RMax_SM2_BF = RMax_V85;
    L_SM2_BF    = L_V85;

    position_SM2_FF[0] = Total_Beam_Length_H2;
    position_SM2_FF[1] = position_of_beam_start[1];
    position_SM2_FF[2] = Total_Beam_Length_V + L_V85/2.0;

    position_SM2_BF[0] = Total_Beam_Length_H2;
    position_SM2_BF[1] = position_of_beam_start[1];
    position_SM2_BF[2] = Total_Beam_Length_V + L_SM2 - L_V85/2.0; 

    rotation_SM2_F[0] =  0.0;
    rotation_SM2_F[1] =  0.0;
    rotation_SM2_F[2] =  0.0;

    // parts-1 
    Rmin_SM2_parts1 = 28.0/2.0;
    Rmax_SM2_parts1 = 34.0/2.0;
    L_SM2_parts1    = L_SM2 - ( L_V85*2.0 );
    
    position_SM2_parts1[0] = Total_Beam_Length_H2;
    position_SM2_parts1[1] = position_of_beam_start[1]; 
    position_SM2_parts1[2] = Total_Beam_Length_V + L_V85 + L_SM2_parts1/2.0;

    // parts-2
    Rmin_SM2_parts2 = 0.0;
    Rmax_SM2_parts2 = 54.0/2.0;
    L_SM2_parts2    = 100.0;
    
    position_SM2_parts2[0] = position_SM2_parts1[0];
    position_SM2_parts2[1] = position_SM2_parts1[1];
    position_SM2_parts2[2] = position_SM2_parts1[2];

    // parts-3
    Rmin_SM2_parts3 = 54.0/2.0;
    Rmax_SM2_parts3 = 60.0/2.0;
    L_SM2_parts3    = 100.0;

    position_SM2_parts3[0] = position_SM2_parts1[0];
    position_SM2_parts3[1] = position_SM2_parts1[1];
    position_SM2_parts3[2] = position_SM2_parts1[2];

    rotation_SM2_comp[0] =  0.0;
    rotation_SM2_comp[1] = 90.0;
    rotation_SM2_comp[2] =  0.0;

    rotation_SM2_duct[0] =  0.0;
    rotation_SM2_duct[1] =  0.0;
    rotation_SM2_duct[2] =  0.0;

    rotation_SM2_body[0] =  0.0;
    rotation_SM2_body[1] = 90.0;
    rotation_SM2_body[2] =  0.0;

    Total_Beam_Length_V += L_SM2;
        
   // No.22 CoreMonitor4 G4Tube     Fragne V85-V85
    RMin_CM4 = 60.4/2.0;
    RMax_CM4 = 66.4/2.0;
    L_CM4    = 200.0 - ( L_V85*2.0 );

    RMin_CM4_FF = RMin_V85;
    RMax_CM4_FF = RMax_V85;
    L_CM4_FF    = L_V85;

    RMin_CM4_BF = RMin_V85;
    RMax_CM4_BF = RMax_V85;
    L_CM4_BF    = L_V85;

    position_CM4[0] = Total_Beam_Length_H2;
    position_CM4[1] = position_of_beam_start[1];
    position_CM4[2] = Total_Beam_Length_V + L_V85+ L_CM4/2.0;

    position_CM4_FF[0] = Total_Beam_Length_H2;
    position_CM4_FF[1] = position_of_beam_start[1];
    position_CM4_FF[2] = Total_Beam_Length_V + L_V85/2.0;

    position_CM4[0] = Total_Beam_Length_H2;
    position_CM4[1] = position_CM4_FF[1];
    position_CM4[2] = Total_Beam_Length_V + L_V85+ L_CM4/2.0;

    position_CM4_BF[0] = Total_Beam_Length_H2;
    position_CM4_BF[1] = position_CM4[1];
    position_CM4_BF[2] = Total_Beam_Length_V + L_V85 + L_CM4 + L_V85/2.0;

    rotation_CM4[0] =  0.0;
    rotation_CM4[1] =  0.0;
    rotation_CM4[2] =  0.0;

    Total_Beam_Length_V += L_CM4 + ( L_V85*2.0 );

   // No.23 Beam Window G4Tube Fragne V85
    RMin_Window = 35.0/2.0;
    RMax_Window = 130.0/2.0;
    L_Window    = 26.0;
    position_Window[0] = Total_Beam_Length_H2;
    position_Window[1] = position_of_beam_start[1];
    position_Window[2] = Total_Beam_Length_V + L_Window/2.0;

    rotation_Window[0] =  0.0;
    rotation_Window[1] =  0.0;
    rotation_Window[2] =  0.0;

    //  No.23.5 TiWindow
    RMin_TiWindow = 0.0;
    RMax_TiWindow = 35.0/2.0;
    L_TiWindow    = 127.0*1E-3; 
    //L_TiWindow    = 200.0*1E-3; // Changed(1/7) for  FC study 2011.04.14
    position_TiWindow[0] = Total_Beam_Length_H2;
    position_TiWindow[1] = position_of_beam_start[1];
    position_TiWindow[2] = Total_Beam_Length_V + 16.0; 
    //position_TiWindow[2] = Total_Beam_Length_V + 26.0; // Changed(2/7) for  FC study 2011.04.14

    rotation_TiWindow[0] =  0.0;
    rotation_TiWindow[1] =  0.0;
    rotation_TiWindow[2] =  0.0;

    L_Window_vac = 16.0 - L_TiWindow/2.0; 
    //L_Window_vac = 26.0 - L_TiWindow/2.0; // Changed(3/7) for  FC study 2011.04.14
    position_Window_vac_z = Total_Beam_Length_V + L_Window_vac/2.0;

    Total_Beam_Length_V += L_Window;  

    cout << " Ti-Window Height (Total_Beam_Length_V : Top of Frange of Beam Exit) = " << Total_Beam_Length_V << endl;

    //  No.23.6 Top Plate
    // Top Plate was removed in 2012.March ...
    size_TopPlate[0] = 900.0; // unit=mm
    size_TopPlate[1] = 240.0; // unit=mm
    size_TopPlate[2] = 25.0;  // unit=mm 
  
    position_TopPlate[0] = Total_Beam_Length_H2 - 20.0;
    position_TopPlate[1] = position_of_beam_start[1];
    position_TopPlate[2] = Total_Beam_Length_V + 10.0 + size_TopPlate[2]/2.0;

    for( int i=0; i<3; i++ ) rotation_TopPlate[i] = 0.0;
   
    RMin_TopPlate_hole = 0.0;
    RMax_TopPlate_hole = 150.0/2.0;
    L_TopPlate_hole    = size_TopPlate[2];

    position_TopPlate_hole[0] = Total_Beam_Length_H2;
    position_TopPlate_hole[1] = position_TopPlate[1]; 
    position_TopPlate_hole[2] = position_TopPlate[2];
    for( int i=0; i<3; i++ ) rotation_TopPlate_hole[i] = 0.0;

    /*  Recover the alignment parameter
        alignment of vertical beam line                                                                                         
     [0] ... shift along x-axis                                                                                              
     [1] ... shift along y-axis */
    Total_Beam_Length_H2-=alignment_vertical_beamline_shift[0];
    position_of_beam_start[1]-=alignment_vertical_beamline_shift[1];

    // No.24 Faraday Cup    
    distance_fc = 52.0; // 10.0; modified by T.Shibata : 2010.02.01
                         /* 52.0mm is real distance 
                            This distance is from frange of Ti-window to Faraday Cup surface.
                            The length from Ti-Window to Faraday cup carbon is 82 mm   
			 */
    first_position_v = Total_Beam_Length_V + distance_fc;
    //first_position_v = Total_Beam_Length_V; // Changed(4/7) for  FC study 2011.04.14

    /*  
     Faraday Cup Alignment ...
    */
     double Total_Beam_Length_H2_FC=Total_Beam_Length_H2+alignment_position_Faradaycup[0];    
     double position_of_beam_start_FC=position_of_beam_start[1]+alignment_position_Faradaycup[1];
     first_position_v+=alignment_position_Faradaycup[2];

    // Al Cylinder
    // parts 1-3
    RMin_FDAlCylinder[0] = 60.0/2.0;
    RMax_FDAlCylinder[0] = 170.0/2.0;
    L_FDAlCylinder[0] = 5.0;

    RMin_FDAlCylinder[1] = 140.0/2.0;
    RMax_FDAlCylinder[1] = 170.0/2.0;
    L_FDAlCylinder[1] = 270.0;
    
    RMin_FDAlCylinder[2] = 0.0;
    RMax_FDAlCylinder[2] = 170.0/2.0;
    L_FDAlCylinder[2] = 5.0;

    for( int i=0; i<3; i++ ) position_FDAlCylinder[i][0] = Total_Beam_Length_H2_FC; //Total_Beam_Length_H2;
    for( int i=0; i<3; i++ ) position_FDAlCylinder[i][1] = position_of_beam_start_FC; //position_of_beam_start[1]; 

    position_FDAlCylinder[0][2] = first_position_v + L_FDAlCylinder[0]/2.0; 
    position_FDAlCylinder[1][2] = first_position_v + L_FDAlCylinder[0] 
                                   + L_FDAlCylinder[1]/2.0;
    position_FDAlCylinder[2][2] = first_position_v + L_FDAlCylinder[0] + L_FDAlCylinder[1]
                                   + L_FDAlCylinder[2]/2.0;

    for( int i=0; i<3; i++ ){
      for( int j=0; j<3; j++ ) rotation_FDAlCylinder[i][j]=0.0; }

    // FD Body
    // parts 1-3 ( Pb ), parts 4 ( C )
    RMin_FDBody[0] = 50.0/2.0;
    RMax_FDBody[0] = 120.0/2.0;
    L_FDBody[0] = 10.0;
    
    RMin_FDBody[1] = 60.0/2.0;
    RMax_FDBody[1] = RMax_FDBody[0];
    L_FDBody[1] = 140.0;

    RMin_FDBody[2] = 0.0;
    RMax_FDBody[2] = RMax_FDBody[0];
    L_FDBody[2] = 100.0;

    RMin_FDBody[3] = 0.0;
    RMax_FDBody[3] = RMin_FDBody[1];
    L_FDBody[3] = L_FDBody[1];
   
    for( int i=0; i<4; i++ ) position_FDBody[i][0] = Total_Beam_Length_H2_FC; //Total_Beam_Length_H2;
    for( int i=0; i<4; i++ ) position_FDBody[i][1] = position_of_beam_start_FC; //position_of_beam_start[1]; 

    delta_fc = 5.0;    
    
    position_FDBody[0][2] = first_position_v + L_FDAlCylinder[0] + delta_fc
                              + L_FDBody[0]/2.0;
    
    //position_FDBody[0][2] = first_position_v + L_FDBody[0]/2.0;  // Changed(5/7) for  FC study 2011.04.14 
    /* Change(6/7) and (7/7) are in DetectorConstruction.cc */

    position_FDBody[1][2] = position_FDBody[0][2] + L_FDBody[0]/2.0
                               + L_FDBody[1]/2.0;
    position_FDBody[2][2] = position_FDBody[1][2] + L_FDBody[1]/2.0
                               + L_FDBody[2]/2.0;
    position_FDBody[3][2] = position_FDBody[1][2];

    for( int i=0; i<3; i++ ){
      for( int j=0; j<3; j++ ) rotation_FDBody[i][j]=0.0; }


    //=============================================
    // Faraday Cup of Copper Body ( =FC2/3/4.... )
    //=============================================
    // FC2-> FC4 in 2012.05.01                                                                                    
    // Front : Carbon                                                                                             
    //       --> Front1/2 thin shield ( 1 and 2 ) Cu or Al or Any other material ?  in 2012.05.01           
    //
    // Creat Faraday Cup 4 in 2012.10.02 
    // Geometry was fixed except position

    // Ground Layer (G4Polycone)
    NumZPlanesGNDLayerFC4 = 8;    

    IrGNDLayerFC4[0] = IrGNDLayerFC4[1] = 0.0; 
    for( int i=2; i<NumZPlanesGNDLayerFC4; i++ ) IrGNDLayerFC4[i] = 90.0/2.0;

    OrGNDLayerFC4[0]=100.0/2.0;
    OrGNDLayerFC4[1]=100.0/2.0;
    OrGNDLayerFC4[2]=100.0/2.0;
    OrGNDLayerFC4[3]=100.0/2.0;
    OrGNDLayerFC4[4]=94.0/2.0;
    OrGNDLayerFC4[5]=94.0/2.0;
    OrGNDLayerFC4[6]=100.0/2.0;
    OrGNDLayerFC4[7]=100.0/2.0;
      
    ZGNDLayerFC4[0] = 0.0;
    ZGNDLayerFC4[1] = 0.1;
    ZGNDLayerFC4[2] = 0.1;
    ZGNDLayerFC4[3] = 10.1;
    ZGNDLayerFC4[4] = 10.1;
    ZGNDLayerFC4[5] = 88.1;
    ZGNDLayerFC4[6] = 88.1;
    ZGNDLayerFC4[7] = 98.1;
    
    position_GNDLayerFC4[0] = Total_Beam_Length_H2 + plist->Positionz_alongbeamline_Faradaycup4();
    position_GNDLayerFC4[1] = position_of_beam_start[1]
                               + plist->Positionz_transverse_Faradaycup4();
    position_GNDLayerFC4[2] = Total_Beam_Length_V
                               + plist->Positionz_Faradaycup4_font1();

    cout << " Position of GNDLayerFC4 = " 
	 << position_GNDLayerFC4[0] << " " 
         << position_GNDLayerFC4[1] << " "
	 << position_GNDLayerFC4[2] << " "
         << endl;  

    for( int i=0; i<3; i++ ){ rotation_GNDLayerFC4[i] = 0.; }

    // Shield Layer (G4Polycone)     
    NumZPlanesShieldLayerFC4 = 14.0;

    IrShieldLayerFC4[0] = 0.0;
    IrShieldLayerFC4[1] = 0.0;
    IrShieldLayerFC4[2] = 75.0/2.0;
    IrShieldLayerFC4[3] = 75.0/2.0;
    IrShieldLayerFC4[4] = 75.0/2.0; 
    IrShieldLayerFC4[5] = 75.0/2.0;
    IrShieldLayerFC4[6] = 20.0/2.0;
    IrShieldLayerFC4[7] = 20.0/2.0;
    IrShieldLayerFC4[8] = 75.0/2.0;
    IrShieldLayerFC4[9] = 75.0/2.0;
    IrShieldLayerFC4[10] = 75.0/2.0;
    IrShieldLayerFC4[11] = 75.0/2.0;
    IrShieldLayerFC4[12] = 15.0/2.0;
    IrShieldLayerFC4[13] = 15.0/2.0;

    OrShieldLayerFC4[0] = 85.0/2.0;
    OrShieldLayerFC4[1] = 85.0/2.0;
    OrShieldLayerFC4[2] = 85.0/2.0;
    OrShieldLayerFC4[3] = 85.0/2.0;
    OrShieldLayerFC4[4] = 79.0/2.0;
    OrShieldLayerFC4[5] = 79.0/2.0;
    OrShieldLayerFC4[6] = 79.0/2.0;    
    OrShieldLayerFC4[7] = 79.0/2.0;    
    OrShieldLayerFC4[8] = 79.0/2.0;    
    OrShieldLayerFC4[9] = 79.0/2.0;    
    OrShieldLayerFC4[10] = 85.0/2.0;    
    OrShieldLayerFC4[11] = 85.0/2.0;    
    OrShieldLayerFC4[12] = 85.0/2.0;    
    OrShieldLayerFC4[13] = 85.0/2.0;    

    ZShieldLayerFC4[0] = 0.0;
    ZShieldLayerFC4[1] = 0.1;
    ZShieldLayerFC4[2] = 0.1;
    ZShieldLayerFC4[3] = 10.1;
    ZShieldLayerFC4[4] = 10.1;
    ZShieldLayerFC4[5] = 64.1;
    ZShieldLayerFC4[6] = 64.1;
    ZShieldLayerFC4[7] = 67.1;
    ZShieldLayerFC4[8] = 67.1;
    ZShieldLayerFC4[9] = 70.1;
    ZShieldLayerFC4[10] = 70.1;
    ZShieldLayerFC4[11] = 80.1;
    ZShieldLayerFC4[12] = 80.1;
    ZShieldLayerFC4[13] = 83.1;

    position_ShieldLayerFC4[0] = position_GNDLayerFC4[0];
    position_ShieldLayerFC4[1] = position_GNDLayerFC4[1]; 
    position_ShieldLayerFC4[2] = position_GNDLayerFC4[2] + 5.0; 
    for( int i=0; i<3; i++ ){ rotation_ShieldLayerFC4[i] = 0.; }

    // Body : Copper : G4Polycone
    NumZPlanesBodyFC4=6;
     
    for( int i=0; i<NumZPlanesBodyFC4; i++ ){ IrBodyFC4[i] = 0.0; }
    OrBodyFC4[0] = 70.0/2.0;
    OrBodyFC4[1] = 70.0/2.0;
    OrBodyFC4[2] = 8.0/2.0;
    OrBodyFC4[3] = 8.0/2.0;
    OrBodyFC4[4] = 20.0/2.0;
    OrBodyFC4[5] = 20.0/2.0;

    ZBodyFC4[0] = 0.0;
    ZBodyFC4[1] = 60.0;
    ZBodyFC4[2] = 60.0;
    ZBodyFC4[3] = 68.0;
    ZBodyFC4[4] = 68.0;
    ZBodyFC4[5] = 71.0;
   
    //RMin_FC4Body = 0.;
    //RMax_FC4Body = plist->R_Faradaycup4_copper();
    //L_FC4Body    = plist->Length_Faradaycup4_copper();

    position_BodyFC4[0] = position_ShieldLayerFC4[0];
    position_BodyFC4[1] = position_ShieldLayerFC4[1];
    position_BodyFC4[2] = position_ShieldLayerFC4[2] + 2.1;
    for( int i=0; i<3; i++ ){ rotation_BodyFC4[i] = 0.; } 
   
    // Isolator-1 : Macor 
    RMin_Iso1FC4=20.0/2.0;                                                                             
    RMax_Iso1FC4=72.0/2.0;                                                                             
    L_Iso1FC4=3.0;
   
    position_Iso1FC4[0] = position_ShieldLayerFC4[0];
    position_Iso1FC4[1] = position_ShieldLayerFC4[1];
    position_Iso1FC4[2] = position_ShieldLayerFC4[2] 
                            + ZShieldLayerFC4[7]
                              + L_Iso1FC4/2.0;
    for( int i=0; i<3; i++ ){ rotation_Iso1FC4[i] = 0.; }
         
    // Isolator-2 : Macor  
    RMin_Iso2FC4=30.0/2.0;
    RMax_Iso2FC4=85.0/2.0;
    L_Iso2FC4=10.0;

    position_Iso2FC4[0] = position_GNDLayerFC4[0];
    position_Iso2FC4[1] = position_GNDLayerFC4[1];
    position_Iso2FC4[2] = position_GNDLayerFC4[2]
                            + ZGNDLayerFC4[5]                            
                              + L_Iso2FC4/2.0;
    for( int i=0; i<3; i++ ){ rotation_Iso2FC4[i] = 0.; }
    
    // Top Plate : A5052
    size_of_TopPlateFC4[0] = 100.0; // mm                     
    size_of_TopPlateFC4[1] = 100.0; // mm  
    size_of_TopPlateFC4[2] = 5.0;   // mm                         

    position_TopPlateFC4[0] = position_GNDLayerFC4[0];
    position_TopPlateFC4[1] = position_GNDLayerFC4[1];
    position_TopPlateFC4[2] = position_GNDLayerFC4[2]
                                + ZGNDLayerFC4[7]
                                  + size_of_TopPlateFC4[2]/2.0;

    for( int i=0; i<3; i++ ){ rotation_TopPlateFC4[i]=0.; }

    //=========================================
    // Faraday Cup 5 
    // made by B.K.Shin-san, in 2013.04.23
    //=========================================
    //Faradaycup 5
    //FC5 No05 body : Copper :G4Polycon

    NumZPlanesBodyFC5=4;
    IrBodyFC5[0]=0.;
    IrBodyFC5[1]=0.;
    IrBodyFC5[2]=4.0/2;;  // M4 --> M5. 
    IrBodyFC5[3]=4.0/2.;  // M4 --> M5. 

    OrBodyFC5[0]=60./2.;
    OrBodyFC5[1]=60./2.;
    OrBodyFC5[2]=60./2.;
    OrBodyFC5[3]=60./2.;

    ZBodyFC5[0]=0;
    ZBodyFC5[1]=48;  // any value is O.K.
    ZBodyFC5[2]=48;  // any value is O.K.
    ZBodyFC5[3]=60;

    position_BodyFC5[0] = Total_Beam_Length_H2      + plist->Faradaycup_shift(0);
    position_BodyFC5[1] = position_of_beam_start[1] + plist->Faradaycup_shift(1);
    position_BodyFC5[2] = Total_Beam_Length_V + plist->Faradaycup5_height() + 17.2; // correctd in '13.10.01
    //position_BodyFC5[2] = position_TiWindow[2] + plist->Faradaycup5_height() + 13.; 
    //position_BodyFC5[2] = position_TiWindow[2] + plist->Faradaycup5_height() + 17.2; // 13 --> 17.2 

    for( int i=0; i<3; i++ ){ rotation_BodyFC5[i] = 0.; }
    
    //FC5 No04 Barrel Suppoter  : Copper :G4Polycon
    NumZPlanesBSFC5=8;
    ZBSFC5[0]=0.;
    ZBSFC5[1]=3.;
    ZBSFC5[2]=3.;
    ZBSFC5[3]=6.;
    ZBSFC5[4]=6.;
    ZBSFC5[5]=67.-2.5-3.;
    ZBSFC5[6]=67.-2.5-3.;
    ZBSFC5[7]=67-2.5;

    IrBSFC5[0]=54.0/2.;
    IrBSFC5[1]=54.0/2.;
    IrBSFC5[2]=60.0/2.;
    IrBSFC5[3]=60.0/2.;
    IrBSFC5[4]=64.5/2.;
    IrBSFC5[5]=64.5/2.;
    IrBSFC5[6]=60./2.;
    IrBSFC5[7]=60./2.;

    OrBSFC5[0]=85./2.;
    OrBSFC5[1]=85./2.;
    OrBSFC5[2]=85./2.;
    OrBSFC5[3]=85./2.;
    OrBSFC5[4]=85./2.;
    OrBSFC5[5]=85./2.;
    OrBSFC5[6]=85./2.;
    OrBSFC5[7]=85./2.;

    position_BSFC5[0] = position_BodyFC5[0];
    position_BSFC5[1] = position_BodyFC5[1];
    position_BSFC5[2] = position_BodyFC5[2]-3;
    for( int i=0; i<3; i++ ){ rotation_BSFC5[i] = 0.; }

    //FC5 No04-Vac1   Barrel Suppoter Vaccum :G4Tubs
    rMaxBSVFC5 = 64.5/2.;
    rMinBSVFC5 = 60./2.;
    HBSVFC5 = (67-2.5-3.-6.)/2.;

    position_BSVFC5[0] = position_BodyFC5[0];
    position_BSVFC5[1] = position_BodyFC5[1];
    position_BSVFC5[2] = position_BodyFC5[2] + 30 + 0.5;

    for( int i=0; i<3; i++ ){ rotation_BSVFC5[i] = 0.; }

    //FC5 No04-Vac2 bottom  Barrel Suppoter Vaccum  :G4Tubs
    rMaxBSVBFC5 = 54/2.;
    rMinBSVBFC5 = 0/2.;
    HBSVBFC5 = 3./2.;

    position_BSVBFC5[0] = position_BodyFC5[0];
    position_BSVBFC5[1] = position_BodyFC5[1];
    position_BSVBFC5[2] = position_BodyFC5[2]-3./2.;

    for( int i=0; i<3; i++ ){ rotation_BSVBFC5[i] = 0.; }

    //FC5 No07 Titan shield (G4Polycon)
    NumZPlanesTSFC5=12;
    ZTSFC5[0]=0.;
    ZTSFC5[1]=0.15; // add componet for bug fixed
    ZTSFC5[2]=0.15;
    ZTSFC5[3]=3.;
    ZTSFC5[4]=3.;
    ZTSFC5[5]=70.-2;
    ZTSFC5[6]=70.-2;
    ZTSFC5[7]=70.-1;
    ZTSFC5[8]=70.-1;
    ZTSFC5[9]=70.;
    ZTSFC5[10]=70.;
    ZTSFC5[11]=72.;

    IrTSFC5[0]=.0/2.;
    IrTSFC5[1]=.0/2.; // add componet for bug fixed
    IrTSFC5[2]=64.0/2.;
    IrTSFC5[3]=65.0/2.;
    IrTSFC5[4]=85.0/2.;
    IrTSFC5[5]=85.5/2.;
    IrTSFC5[6]=60./2.;
    IrTSFC5[7]=60./2.;
    IrTSFC5[8]=61./2.;
    IrTSFC5[9]=61./2.;
    IrTSFC5[10]=85./2.;
    IrTSFC5[11]=85./2.;

    OrTSFC5[0]=90./2.;
    OrTSFC5[1]=90./2.; // add componet for bug fixed 
    OrTSFC5[2]=90./2.;
    OrTSFC5[3]=90./2.;
    OrTSFC5[4]=90./2.;
    OrTSFC5[5]=90./2.;
    OrTSFC5[6]=90./2.;
    OrTSFC5[7]=90./2.;
    OrTSFC5[8]=90./2.;
    OrTSFC5[9]=90./2.;
    OrTSFC5[10]=90./2.;
    OrTSFC5[11]=90./2.;

    position_TSFC5[0] = position_BodyFC5[0];
    position_TSFC5[1] = position_BodyFC5[1];
    position_TSFC5[2] = position_BodyFC5[2]-6;

    for( int i=0; i<3; i++ ){ rotation_TSFC5[i] = 0.; }

    //FC5 No07-Vac    Titan sheild Vaccum :G4Cons
    brMaxTSVFC5 = 64.0/2.;
    brMinTSVFC5 = 0./2.;
    trMaxTSVFC5 = 65.0/2.;
    trMinTSVFC5 = 0./2.;
    HTSVFC5 = (3.-0.15)/2.;
    position_TSVFC5[0] = position_BodyFC5[0];
    position_TSVFC5[1] = position_BodyFC5[1];
    position_TSVFC5[2] = position_BodyFC5[2]-3-HTSVFC5;

    for( int i=0; i<3; i++ ){ rotation_TSVFC5[i] = 0.; }

    //FC5 No06 Bottom Supporter (BotS) (G4Polycon)
    NumZPlanesBotSFC5=4;
    ZBotSFC5[0]=0.;
    ZBotSFC5[1]=4.;
    ZBotSFC5[2]=4.;
    ZBotSFC5[3]=9.;

    IrBotSFC5[0]=68./2.;
    IrBotSFC5[1]=68./2.;
    IrBotSFC5[2]=90./2.;
    IrBotSFC5[3]=90./2.;

    OrBotSFC5[0]=95./2.;
    OrBotSFC5[1]=95./2.;
    OrBotSFC5[2]=95./2.;
    OrBotSFC5[3]=95./2.;

    position_BotSFC5[0] = position_BodyFC5[0];
    position_BotSFC5[1] = position_BodyFC5[1];
    position_BotSFC5[2] = position_BodyFC5[2]-10;
    for( int i=0; i<3; i++ ){ rotation_BotSFC5[i] = 0.; }

    //FC5 No15 Top Supporter (TopS) (G4Polycon)
    NumZPlanesTopSFC5=6;
    ZTopSFC5[0]=0.;
    ZTopSFC5[1]=3.;
    ZTopSFC5[2]=3.;
    ZTopSFC5[3]=22.;
    ZTopSFC5[4]=22.;
    ZTopSFC5[5]=25.;

    IrTopSFC5[0]=50./2.;
    IrTopSFC5[1]=50./2.;
    IrTopSFC5[2]=50./2.;
    IrTopSFC5[3]=50./2.;
    IrTopSFC5[4]=50./2.;
    IrTopSFC5[5]=50./2.;

    OrTopSFC5[0]=59./2.;
    OrTopSFC5[1]=59./2.;
    OrTopSFC5[2]=61./2.;
    OrTopSFC5[3]=61./2.;
    OrTopSFC5[4]=59./2.;
    OrTopSFC5[5]=59./2.;

    position_TopSFC5[0] = position_BodyFC5[0];
    position_TopSFC5[1] = position_BodyFC5[1];
    position_TopSFC5[2] = position_BodyFC5[2]+60;

    for( int i=0; i<3; i++ ){ rotation_TopSFC5[i] = 0.; }

    //FC5 No15-Vac  Top  Suppoter Vaccum  :G4Tubs
    rMaxTopSVFC5 = 50./2.;
    rMinTopSVFC5 = 4./2.;
    HTopSVFC5 = 25./2.;

    position_TopSVFC5[0] = position_BodyFC5[0];
    position_TopSVFC5[1] = position_BodyFC5[1];
    position_TopSVFC5[2] = position_BodyFC5[2]+60+HTopSVFC5;

    for( int i=0; i<3; i++ ){ rotation_TopSVFC5[i] = 0.; }

    //FC5 No02 & 01  CF90  (CF90)   (G4Polycon)
    NumZPlanesCF90FC5=6;
    ZCF90FC5[0]=0.;
    ZCF90FC5[1]=2.;
    ZCF90FC5[2]=2.;
    ZCF90FC5[3]=21.;
    ZCF90FC5[4]=21.;
    ZCF90FC5[5]=35.;

    IrCF90FC5[0]=65./2.;
    IrCF90FC5[1]=61./2.;
    IrCF90FC5[2]=61./2.;
    IrCF90FC5[3]=61./2.;
    IrCF90FC5[4]=4./2.;
    IrCF90FC5[5]=4./2.;

    OrCF90FC5[0]=85./2.;
    OrCF90FC5[1]=85./2.;
    OrCF90FC5[2]=90./2.;
    OrCF90FC5[3]=90./2.;
    OrCF90FC5[4]=90./2.;
    OrCF90FC5[5]=90./2.;

    position_CF90FC5[0] = position_BodyFC5[0];
    position_CF90FC5[1] = position_BodyFC5[1];
    position_CF90FC5[2] = position_BodyFC5[2]+70-6;

    for( int i=0; i<3; i++ ){ rotation_CF90FC5[i] = 0.; }

    //FC5 No02 & 01  bTopS  (bTopS)   (G4Polycon)
    NumZPlanesbTopSFC5=4;
    ZbTopSFC5[0]=0.;
    ZbTopSFC5[1]=10.;
    ZbTopSFC5[2]=10.;
    ZbTopSFC5[3]=15.;

    IrbTopSFC5[0]=90./2.;
    IrbTopSFC5[1]=90./2.;
    IrbTopSFC5[2]=85./2.;
    IrbTopSFC5[3]=85./2.;

    OrbTopSFC5[0]=95./2.;
    OrbTopSFC5[1]=95./2.;
    OrbTopSFC5[2]=95./2.;
    OrbTopSFC5[3]=95./2.;

    position_bTopSFC5[0] = position_BodyFC5[0];
    position_bTopSFC5[1] = position_BodyFC5[1];
    position_bTopSFC5[2] = position_BodyFC5[2]+64.+20.+12.;

    for( int i=0; i<3; i++ ){ rotation_bTopSFC5[i] = 0.; }

    //Copper chamber  No11 (CP)
    NumZPlanesCPFC5=10;
    ZCPFC5[0]=0.;
    ZCPFC5[1]=4.; 
    ZCPFC5[2]=4.; 
    ZCPFC5[3]=4.2;
    ZCPFC5[4]=4.2;
    ZCPFC5[5]=7.2;
    ZCPFC5[6]=7.2;
    ZCPFC5[7]=130+4.2;
    ZCPFC5[8]=130+4.2;
    ZCPFC5[9]=140+4.2;

    IrCPFC5[0]=75./2.;
    IrCPFC5[1]=75./2.;
    IrCPFC5[2]=0./2.;
    IrCPFC5[3]=0./2.;
    IrCPFC5[4]=80./2.;
    IrCPFC5[5]=80./2.;
    IrCPFC5[6]=95./2.;
    IrCPFC5[7]=95./2.;
    IrCPFC5[8]=95./2.;
    IrCPFC5[9]=95./2.;

    OrCPFC5[0]=99./2.;
    OrCPFC5[1]=99./2.;
    OrCPFC5[2]=99./2.;
    OrCPFC5[3]=99./2.;
    OrCPFC5[4]=99./2.;
    OrCPFC5[5]=99./2.;
    OrCPFC5[6]=99./2.;
    OrCPFC5[7]=99./2.;
    OrCPFC5[8]=110./2.;
    OrCPFC5[9]=110./2.;

    position_CPFC5[0] =  position_BodyFC5[0];
    position_CPFC5[1] =  position_BodyFC5[1];
    position_CPFC5[2] =  position_BodyFC5[2] - 17.2;
    for( int i=0; i<3; i++ ){ rotation_CPFC5[i] = 0.; }

    //Parts No 12 Top Plate TP
    XTPFC5=110./2.;
    YTPFC5=195.5/2.;
    ZTPFC5=5./2.;

    position_TPFC5[0] =  position_BodyFC5[0];
    position_TPFC5[1] =  position_BodyFC5[1] +195./2.-56;
    position_TPFC5[2] =  position_BodyFC5[2] -13.+140.+ZTPFC5;
    for( int i=0; i<3; i++ ){ rotation_TPFC5[i] = 0.; }

    //FC5 FeedThrough :G4Tubs
    rMaxFTFC5 = 4./2.;
    rMinFTFC5 = 0/2.;
    HFTFC5 = (32.)/2.;

    position_FTFC5[0] =   position_BodyFC5[0];
    position_FTFC5[1] =  position_BodyFC5[1];
    position_FTFC5[2] =  position_BodyFC5[2]+60.-12.+HFTFC5;
    for( int i=0; i<3; i++ ){ rotation_FTFC5[i] = 0.; }

    //FC5 FeedThrough2 :G4Polycon
    rMaxFT2FC5 = 4./2.;
    rMinFT2FC5 = 0/2.;
    HFT2FC5 = (20.)/2.;

    position_FT2FC5[0] =  position_BodyFC5[0];
    position_FT2FC5[1] =  position_BodyFC5[1];
    position_FT2FC5[2] =   position_FTFC5[2]+HFTFC5+HFT2FC5;
    for( int i=0; i<3; i++ ){ rotation_FT2FC5[i] = 0.; }

    //FC5 VaccumDuct Vacuum :G4Tubs
    rMaxVDVFC5 = 18./2.;
    rMinVDVFC5 = 0/2.;
    HVDVFC5 = 113./2.;

    position_VDVFC5[0] =   position_BodyFC5[0];
    position_VDVFC5[1] =  position_BodyFC5[1]-20;
    position_VDVFC5[2] =  position_BodyFC5[2]+60.+25.+HVDVFC5;
    for( int i=0; i<3; i++ ){ rotation_VDVFC5[i] = 0.; }

    //FC5 VaccumDuct ::G4Tubs
    rMaxVDFC5 = 20./2.;
    rMinVDFC5 = 18./2.;
    HVDFC5 = 100./2.;

    position_VDFC5[0] =   position_BodyFC5[0];
    position_VDFC5[1] =    position_BodyFC5[1]-20;
    position_VDFC5[2] =    position_CF90FC5[2]+35.+HVDFC5;
    for( int i=0; i<3; i++ ){ rotation_VDFC5[i] = 0.; }

    //===========================================
    // Screen Monitor 3
    //=========================================== 
    size_of_ScreenMonitor3[0] = 100.0; // mm         
    size_of_ScreenMonitor3[1] = 100.0; // mm
    size_of_ScreenMonitor3[2] = 0.5;   // mm

    position_ScreenMonitor3[0] = Total_Beam_Length_H2;
    position_ScreenMonitor3[1] = position_of_beam_start[1];
    position_ScreenMonitor3[2] = Total_Beam_Length_V + plist->Positionz_ScreenMonitor3() +
                                 size_of_ScreenMonitor3[2]/2.0;
    
    rotation_ScreenMonitor3[0] = -45.0;
    rotation_ScreenMonitor3[1] = 0.0;
    rotation_ScreenMonitor3[2] = 0.0;
    
    //===========================================
    // Beam Attenuator
    //===========================================
    RMin_BeamAttenuator = 0;
    //RMax_BeamAttenuator = 55.0/2.0; 
    //RMax_BeamAttenuator = 70.0/2.0; 
    RMax_BeamAttenuator = 60.0/2.0; 
    L_BeamAttenuator    = plist->Length_Beam_Attenuator();

    position_BeamAttenuator[0] = Total_Beam_Length_H2;
    position_BeamAttenuator[1] = position_of_beam_start[1];
    position_BeamAttenuator[2] = Total_Beam_Length_V + L_BeamAttenuator/2.0;

    for( int i=0; i<3; i++ ){ rotation_BeamAttenuator[i] = 0.; } 

    //===========================================
    // Pb Collimator added in 2012.09.24
    //===========================================
    RMin_PbCollimator = plist->Dinner_Pb_Collimator()/2.0;
    RMax_PbCollimator = plist->Douter_Pb_Collimator()/2.0;
    L_PbCollimator    = plist->Length_Pb_Collimator();

    position_PbCollimator[0] = Total_Beam_Length_H2;
    position_PbCollimator[1] = position_of_beam_start[1];
    position_PbCollimator[2] = Total_Beam_Length_V + 2.0 + L_PbCollimator/2.0; 

    for( int i=0; i<3; i++ ){ rotation_PbCollimator[i] = 0.; }

    //==================================
    // Container Injection Hole 
    //==================================
    size_of_ELS_injection_hole[0] = 770.0;
    size_of_ELS_injection_hole[1] = 770.0;
    size_of_ELS_injection_hole[2] = 500.0; // not be used any more...

    position_of_ELS_injection_hole[0] = Total_Beam_Length_H2;
    position_of_ELS_injection_hole[1] = position_of_beam_start[1];
    position_of_ELS_injection_hole[2] = position_of_ELS[2] 
                                        + outer_size_of_ELScontainer[2]/2.0;

    //==================================
    // Cover Box
    //==================================
  
    size_of_coverbox_outer[0] = size_of_ELS_injection_hole[0];
    size_of_coverbox_outer[1] = size_of_ELS_injection_hole[1];
    size_of_coverbox_outer[2] = 814.0; // 485( Coverbox height )+20(Thickness of boad)+309(Inner)

    thickness_of_coverbox = 5.0;    
    size_of_coverbox_inner[0] = size_of_coverbox_outer[0] - thickness_of_coverbox*2.0;
    size_of_coverbox_inner[1] = size_of_coverbox_outer[2] - thickness_of_coverbox*2.0;
    size_of_coverbox_inner[2] = size_of_coverbox_outer[2];

    thickness_of_coverbox_board = 20.0;
    size_of_coverbox_board[0] = size_of_coverbox_inner[0];
    size_of_coverbox_board[1] = size_of_coverbox_inner[1];
    size_of_coverbox_board[2] = thickness_of_coverbox_board;

    rmin_coverbox_boardhole = 0.0;
    rmax_coverbox_boardhole = 200.0/2.0;
    l_coverbox_boardhole    = thickness_of_coverbox_board;

    position_of_coverbox_outer[0] = Total_Beam_Length_H2;
    position_of_coverbox_outer[1] = position_of_beam_start[1];
    position_of_coverbox_outer[2] = position_of_ELS[2]  
                                      + outer_size_of_ELScontainer[2]/2.0
                                        - thick_of_ELScontainer_wall_roof
                                           - 30.0
                                            + size_of_coverbox_outer[2]/2.0;

    for( int i=0; i<3; i++ ) position_of_coverbox_inner[i] = position_of_coverbox_outer[i];
    
    zbias_of_coverbox_board = 309.0;
    position_of_coverbox_board[0] = position_of_coverbox_outer[0];
    position_of_coverbox_board[1] = position_of_coverbox_outer[1];
    position_of_coverbox_board[2] = position_of_ELS[2]
                                      + outer_size_of_ELScontainer[2]/2.0
                                        - thick_of_ELScontainer_wall_roof
                                          - 30.0
                                           + zbias_of_coverbox_board
                                             + size_of_coverbox_board[2]/2.0;
 
    position_of_coverbox_boardhole[0] = position_of_coverbox_board[0] 
                                           + alignment_coverbox_boardhole_shift[0];
    position_of_coverbox_boardhole[1] = position_of_coverbox_board[1] 
                                           + alignment_coverbox_boardhole_shift[1];
    position_of_coverbox_boardhole[2] = position_of_coverbox_board[2];

    //==========================================================
    // Virtual Test Chamber added in 2011.12.13
    //   modified in 2013.04.17, to check backscattering effect
    //                     all of gemetory parameters are fixed
    //==========================================================
    // Cylinder( Material ) parts 1-3                  
    // parts-1 : Front Panel
    // parts-2 : Cylinder
    // parts-3 : Rear Panel
    // parts-4 : Inner cylinder -1  --> not be used in 2013.04.17
    // parts-5 : Inner cylinder -2  --> not be used in 2013.04.17

    double TargetD=60.0;  // 10000.0;
    RMin_VirtualChamberCylinder[0] = 0.0;
    RMin_VirtualChamberCylinder[1] = TargetD/2.0 + 1.0 + 0.9;
    RMin_VirtualChamberCylinder[2] = TargetD/2.0 + 1.0;      
    RMin_VirtualChamberCylinder[3] = RMin_VirtualChamberCylinder[2];      
    RMin_VirtualChamberCylinder[4] = 0.0;

    RMax_VirtualChamberCylinder[0] = TargetD/2.0 + 1.0 + 1.0; 
    RMax_VirtualChamberCylinder[1] = RMax_VirtualChamberCylinder[0];    
    RMax_VirtualChamberCylinder[2] = RMax_VirtualChamberCylinder[0];
    RMax_VirtualChamberCylinder[3] = RMax_VirtualChamberCylinder[0];
    RMax_VirtualChamberCylinder[4] = RMax_VirtualChamberCylinder[0];

    L_VirtualChamberCylinder[0] = 0.5;
    L_VirtualChamberCylinder[1] = 1.0;
    L_VirtualChamberCylinder[2] = 2.0 + 10000.0/2.0; 
    L_VirtualChamberCylinder[3] =       10000.0/2.0 + 2.0; 
    L_VirtualChamberCylinder[4] = 1.0;

    for( int i=0; i<3; i++ ) position_VirtualChamberCylinder[i][0] = Total_Beam_Length_H2;
    for( int i=0; i<3; i++ ) position_VirtualChamberCylinder[i][1] = position_of_beam_start[1];

    position_VirtualChamberCylinder[0][2] = plist->Positionz_virtual_chamber() 
                                             + L_VirtualChamberCylinder[0]/2.0;

    position_VirtualChamberCylinder[1][2] = position_VirtualChamberCylinder[0][2]
                                             + L_VirtualChamberCylinder[0]/2.0 
                                               + L_VirtualChamberCylinder[1]/2.0;     

    position_VirtualChamberCylinder[2][2] = position_VirtualChamberCylinder[1][2]
                                             + L_VirtualChamberCylinder[1]/2.0
                                               + L_VirtualChamberCylinder[2]/2.0;

    position_VirtualChamberCylinder[3][2] = position_VirtualChamberCylinder[2][2]
                                             + L_VirtualChamberCylinder[2]/2.0
                                               + L_VirtualChamberCylinder[3]/2.0; 

    position_VirtualChamberCylinder[4][2] = position_VirtualChamberCylinder[3][2]
                                               + L_VirtualChamberCylinder[3]/2.0
                                                //+ plist->Gap_virtual_chamber() 
                                                 + L_VirtualChamberCylinder[4]/2.0;

    for( int i=0; i<5; i++ ){
      for( int j=0; j<3; j++ ){ rotation_VirtualChamberCylinder[i][j]=0.0; } 
    }
    // Cylinder is OK ...
 
    // Cu Target : added in 2013.04.17 
    RMin_VirtualChamberTarget = 0.0; 
    RMax_VirtualChamberTarget = TargetD/2.0; 
    L_VirtualChamberTarget    = 10000.0;  
   
    position_VirtualChamberTarget[0] = position_VirtualChamberCylinder[0][0];
    position_VirtualChamberTarget[1] = position_VirtualChamberCylinder[0][1];
    position_VirtualChamberTarget[2] = ( position_VirtualChamberCylinder[2][2] + 
                                         position_VirtualChamberCylinder[3][2] )/2.0;
    
    rotation_VirtualChamberTarget[0] = 0.0;
    rotation_VirtualChamberTarget[1] = 0.0;
    rotation_VirtualChamberTarget[2] = 0.0;

    // Body : Inner side : Air ( N2 or ... ) parts 1-4
    RMin_VirtualChamberBody[0] = RMax_VirtualChamberCylinder[3];
    RMin_VirtualChamberBody[1] = 0.0;
    RMin_VirtualChamberBody[2] = 0.0; 
    RMin_VirtualChamberBody[3] = 0.0;

    RMax_VirtualChamberBody[0] = plist->R_virtual_chamber_innerhole();
    RMax_VirtualChamberBody[1] = plist->R_virtual_chamber_injectionhole();
    RMax_VirtualChamberBody[2] = RMin_VirtualChamberBody[0];
    RMax_VirtualChamberBody[3] = plist->R_virtual_chamber_injectionhole();
    
    L_VirtualChamberBody[0]    = L_VirtualChamberCylinder[1];
    L_VirtualChamberBody[1]    = L_VirtualChamberCylinder[3];
    L_VirtualChamberBody[2]    = plist->Gap_virtual_chamber();
    L_VirtualChamberBody[3]    = L_VirtualChamberCylinder[3];

    for( int i=0; i<4; i++ ){
      position_VirtualChamberBody[i][0] = Total_Beam_Length_H2;
      position_VirtualChamberBody[i][1] = position_of_beam_start[1];
    }

    position_VirtualChamberBody[0][2] = position_VirtualChamberCylinder[1][2];
    position_VirtualChamberBody[1][2] = position_VirtualChamberCylinder[3][2];
    position_VirtualChamberBody[2][2] = position_VirtualChamberCylinder[1][2];
    position_VirtualChamberBody[3][2] = position_VirtualChamberCylinder[4][2];

    for( int i=0; i<4; i++ ){ 
      for( int j=0; j<3; j++ ){ 
          rotation_VirtualChamberBody[i][j]=0; 
      }
    }

    //==================================
    // Beam injection Position
    //==================================
    beam_injection_position[0] = Total_Beam_Length_H2;
    beam_injection_position[1] = position_of_beam_start[1];
    beam_injection_position[2] = 0.0; 


    //==================================
    // Define ICE
    //==================================
    ice_position[0] = Total_Beam_Length_H2;
    ice_position[1] = position_of_beam_start[1];
    //ice_position[2] = 7500.0; //5000.0; // mm -- this needs to be above concrete walls.  See picture from chat
    ice_position[2] = 10000; //7500.0; //5000.0; // mm -- this needs to be above concrete walls.  See picture from chat

    ice_size[0] = 500.0;   //mm
    ice_size[1] = 4000.0;  //mm
    ice_size[2] = 10000.0; //5000.0;   //mm

    ice_rotation[0] = 0.; // no rotation for now
    ice_rotation[1] = 0.;
    ice_rotation[2] = 0.;
    
    
    
    // -------------------------
    inputBeamLineG4TubeList();
    // -------------------------
    
    // Quadrupole Magnet 1 Field
    // X-axis
    quadrupole_magnet1_field_region[0][0] = position_QM_FF[0]-L_V85/2.0 + 75.75;
    quadrupole_magnet1_field_region[0][1] = position_QM_FF[0]-L_V85/2.0 + 75.75 + 221.0;
    // Y-axis
    quadrupole_magnet1_field_region[1][0] = position_of_beam_start[1] - 50.0;
    quadrupole_magnet1_field_region[1][1] = position_of_beam_start[1] + 50.0;
    // Z-axis
    quadrupole_magnet1_field_region[2][0] = position_of_beam_start[2] - 50.0;
    quadrupole_magnet1_field_region[2][1] = position_of_beam_start[2] + 50.0;

    // Quadrupole Magnet 2 Field
    // X-axis
    quadrupole_magnet2_field_region[0][0] = position_QM_FF[0]-L_V85/2.0 + 75.75 + 221.0 + 79.0;
    quadrupole_magnet2_field_region[0][1] = position_QM_FF[0]-L_V85/2.0 + 75.75 + 221.0 + 79.0 + 221.0;
    // Y-axis
    quadrupole_magnet2_field_region[1][0] = position_of_beam_start[1] - 50.0;
    quadrupole_magnet2_field_region[1][1] = position_of_beam_start[1] + 50.0;
    // Z-axis
    quadrupole_magnet2_field_region[2][0] = position_of_beam_start[2] - 50.0;
    quadrupole_magnet2_field_region[2][1] = position_of_beam_start[2] + 50.0;

    // Bending Magnetic Field
    // X-axis
    bending_magnet_field_region[0][0] = Total_Beam_Length_H2 - 200.0;
    bending_magnet_field_region[0][1] = Total_Beam_Length_H2 + 50.0;  
    
    // Y-axis
    bending_magnet_field_region[1][0] = position_of_beam_start[1] - 12.0;
    bending_magnet_field_region[1][1] = position_of_beam_start[1] + 12.0;

    // Z-axis
    bending_magnet_field_region[2][0] = position_of_beam_start[2] - 50.0;
    bending_magnet_field_region[2][1] = position_of_beam_start[2] + 200.0;
 
    /* alignment of BM York & Coil
      0 : shift along x-axis
      1 : shift along y-axis
      2 : shift along z-axis
    */
    for( int i=0; i<3; i++ ){
      for( int j=0; j<2; j++ ){
           bending_magnet_field_region[i][j]+=alignment_position_BMYork_Coil[i];
      }
    }

}
//------------------------------------------------
void ELSParameters::shielding_parameter_initialize(void)
{

  //---------------------------------------
  // Pb Block
  //---------------------------------------

  // PB Block tower
  //  0: Forward
  //  1: Backward
  //  2: Side
  //  3: Side
  //  5: Straight Line

  standard_pb_L=200.0;
  standard_pb_W=100.0;
  standard_pb_H=50.0;    

  additional_pb_hight = 300.0;

  // 0: Forward
  size_of_PbTower[0][0] = standard_pb_H;     //=50.0
  size_of_PbTower[0][1] = standard_pb_W*8.0; //=800.0
  //  size_of_PbTower[0][2] = standard_pb_L*5.0+additional_pb_hight; //=1300.0
  size_of_PbTower[0][2] = standard_pb_L*5.0; // Additional PB Block is not need.

  position_of_PbTower[0][0] = Total_Beam_Length_H2 
                              - ( 500.0 - size_of_PbTower[0][0]/2.0 ) ;
  position_of_PbTower[0][1] = position_of_beam_start[1];
  position_of_PbTower[0][2] = position_of_beam_start[2]
                               -(1173-994.3) + standard_pb_L*2.0 
                                + size_of_PbTower[0][2]/2.0;

  // 1: Backward
  for( int i=0; i<3; i++ ) size_of_PbTower[1][i] =  size_of_PbTower[0][i];
  position_of_PbTower[1][0] = Total_Beam_Length_H2
                              + ( 500.0 - size_of_PbTower[1][0]/2.0 ) ;
  position_of_PbTower[1][1] = position_of_PbTower[0][1];
  position_of_PbTower[1][2] = position_of_PbTower[0][2];

  // 2: Side 
  size_of_PbTower[2][0] = standard_pb_W*10.0;
  size_of_PbTower[2][1] = standard_pb_H;
  //size_of_PbTower[2][2] = standard_pb_L*7.0 + additional_pb_hight; //=1700.0
  size_of_PbTower[2][2] = standard_pb_L*7.0; // Additional PB Block is not need.

  position_of_PbTower[2][0] = Total_Beam_Length_H2;
  position_of_PbTower[2][1] = position_of_beam_start[1] 
                              - ( 400.0 + standard_pb_H/2.0 );
  position_of_PbTower[2][2] = position_of_beam_start[2]
                              -(1173-994.3) + size_of_PbTower[2][2]/2.0;

  // 3: Side
  for( int i=0; i<3; i++ ) size_of_PbTower[3][i] =  size_of_PbTower[2][i];
  position_of_PbTower[3][0] = position_of_PbTower[2][0];
  position_of_PbTower[3][1] = position_of_beam_start[1]
                              + ( 400.0 + standard_pb_H/2.0 );
  position_of_PbTower[3][2] = position_of_PbTower[2][2];

  // 4: Straight Line
  size_of_PbTower[4][0] = standard_pb_L;
  size_of_PbTower[4][1] = standard_pb_W*8.0; // = 800.0
  size_of_PbTower[4][2] = standard_pb_H*14.0; // = 700.0
  position_of_PbTower[4][0] = Total_Beam_Length_H2 
                              + 479.0 + size_of_PbTower[4][0]/2.0; 
  position_of_PbTower[4][1] = position_of_beam_start[1];
  position_of_PbTower[4][2] = position_of_beam_start[2]
                              - (1173-994.3) + size_of_PbTower[4][2]/2.0;


  // BM Duct Pb Block
  //  0: Forward
  //  1: Backward
  
  // 0: Forward
  // parts-1
  Rmin_PbBlock_BM_parts1[0] = 115.0;
  Rmax_PbBlock_BM_parts1[0] = 165.0;
  L_PbBlock_BM_parts1[0] = 22.0;
  Phi_start_PbBlock_BM_parts1[0] = 270.0;
  Phi_delta_PbBlock_BM_parts1[0] = 90.0;

  position_of_PbBlock_BM_parts1[0][0] = Total_Beam_Length_H2 
                                        - 35.0 - Rmax_PbBlock_BM_parts1[0];
  position_of_PbBlock_BM_parts1[0][1] = position_of_beam_start[1];
  position_of_PbBlock_BM_parts1[0][2] = position_of_beam_start[2] 
                                        + 35.0 + Rmax_PbBlock_BM_parts1[0];
  rotation_of_PbBlock_BM_parts1[0][0] = 90.0;
  rotation_of_PbBlock_BM_parts1[0][1] = 0.0;
  rotation_of_PbBlock_BM_parts1[0][2] = 0.0;
  // parts-2
  size_of_PbBlock_BM_parts2[0][0] = 50.0;
  size_of_PbBlock_BM_parts2[0][1] = 22.0;
  size_of_PbBlock_BM_parts2[0][2] = 120.0;
  position_of_PbBlock_BM_parts2[0][0] = Total_Beam_Length_H2
                                     - 35.0 - size_of_PbBlock_BM_parts2[0][0]/2.0;
  position_of_PbBlock_BM_parts2[0][1] = position_of_beam_start[1];
  position_of_PbBlock_BM_parts2[0][2] = position_of_beam_start[2]
                                      + 35.0 + Rmax_PbBlock_BM_parts1[0]
                                      + size_of_PbBlock_BM_parts2[0][2]/2.0;
  // parts-3
  size_of_PbBlock_BM_parts3[0][0] = 120.0;
  size_of_PbBlock_BM_parts3[0][1] = 22.0;
  size_of_PbBlock_BM_parts3[0][2] = 50.0;
  position_of_PbBlock_BM_parts3[0][0] = Total_Beam_Length_H2
                                        - 35.0 - Rmax_PbBlock_BM_parts1[0]
                                        - size_of_PbBlock_BM_parts3[0][0]/2.0;
  position_of_PbBlock_BM_parts3[0][1] = position_of_beam_start[1];
  position_of_PbBlock_BM_parts3[0][2] = position_of_beam_start[2]
                                        + 35.0 + size_of_PbBlock_BM_parts3[0][2]/2.0;
   
  // 1: Backward
  // parts-1
  Rmin_PbBlock_BM_parts1[1] = 0.0;
  Rmax_PbBlock_BM_parts1[1] = 50.0;
  L_PbBlock_BM_parts1[1] = 22.0;
  Phi_start_PbBlock_BM_parts1[1] = 180.0;
  Phi_delta_PbBlock_BM_parts1[1] = 90.0;

  position_of_PbBlock_BM_parts1[1][0] = Total_Beam_Length_H2
                                        + 35.0 + Rmax_PbBlock_BM_parts1[1];
  position_of_PbBlock_BM_parts1[1][1] = position_of_beam_start[1];
  position_of_PbBlock_BM_parts1[1][2] = position_of_beam_start[2]
                                        + 35.0 + Rmax_PbBlock_BM_parts1[1];
  rotation_of_PbBlock_BM_parts1[1][0] = 90.0;
  rotation_of_PbBlock_BM_parts1[1][1] = 0.0;
  rotation_of_PbBlock_BM_parts1[1][2] = 0.0;
  // parts-2
  size_of_PbBlock_BM_parts2[1][0] = 50.0;
  size_of_PbBlock_BM_parts2[1][1] = 22.0;
  size_of_PbBlock_BM_parts2[1][2] = 235.0;
  position_of_PbBlock_BM_parts2[1][0] = Total_Beam_Length_H2
                                         + 35.0 + size_of_PbBlock_BM_parts2[1][0]/2.0;
  position_of_PbBlock_BM_parts2[1][1] = position_of_beam_start[1];
  position_of_PbBlock_BM_parts2[1][2] = position_of_beam_start[2]
                                      + 35.0 + Rmax_PbBlock_BM_parts1[1]
                                      + size_of_PbBlock_BM_parts2[1][2]/2.0;
  // parts-3
  size_of_PbBlock_BM_parts3[1][0] = 85.0;
  size_of_PbBlock_BM_parts3[1][1] = 22.0;
  size_of_PbBlock_BM_parts3[1][2] = 50.0;
  position_of_PbBlock_BM_parts3[1][0] = Total_Beam_Length_H2
                                        + 35.0 + Rmax_PbBlock_BM_parts1[1]
                                        + size_of_PbBlock_BM_parts3[1][0]/2.0;
  position_of_PbBlock_BM_parts3[1][1] = position_of_beam_start[1];
  position_of_PbBlock_BM_parts3[1][2] = position_of_beam_start[2]
                                        + 35.0 + size_of_PbBlock_BM_parts3[1][2]/2.0;


  //---------------------------------------
  // Concrete Block
  //---------------------------------------
  numConcrete = 9;

  thick_concrete = 2.0*ft2mm;
  height_concrete = 12.0*ft2mm;
  posiz_concrete  = position_of_concretepad[2] 
                     + size_of_concretepad[2]/2.0 + height_concrete/2.0; 

  for( int i=0; i < numConcrete ; i++ ){
    size_of_ConcreteBlock[i][2] = height_concrete;
    position_of_ConcreteBlock[i][2] = posiz_concrete; 
  }

  // block #1
  size_of_ConcreteBlock[0][0] = thick_concrete;
  size_of_ConcreteBlock[0][1] = 10.0*ft2mm;
  position_of_ConcreteBlock[0][0] = -158.0 - 100.0
                                    - size_of_ConcreteBlock[0][0]/2.0; 
  position_of_ConcreteBlock[0][1] = 158.0 + 200.0 
                                    - size_of_ConcreteBlock[0][1]/2.0;
  // block #2
  size_of_ConcreteBlock[1][0] = 6.0*ft2mm;
  size_of_ConcreteBlock[1][1] = thick_concrete;
  position_of_ConcreteBlock[1][0] = position_of_ConcreteBlock[0][0] 
                                     - size_of_ConcreteBlock[0][0]/2.0 
                                     + size_of_ConcreteBlock[1][0]/2.0;
  position_of_ConcreteBlock[1][1] = position_of_ConcreteBlock[0][1]
                                     - size_of_ConcreteBlock[0][1]/2.0 
                                     - size_of_ConcreteBlock[1][1]/2.0; 
  // block #3
  size_of_ConcreteBlock[2][0] = thick_concrete;
  size_of_ConcreteBlock[2][1] = 4.0*ft2mm;
  position_of_ConcreteBlock[2][0] = position_of_ConcreteBlock[1][0] 
                                     + size_of_ConcreteBlock[1][0]/2.0
                                     - size_of_ConcreteBlock[2][0]/2.0;
  position_of_ConcreteBlock[2][1] = position_of_ConcreteBlock[1][1] 
                                     - size_of_ConcreteBlock[1][1]/2.0
                                     - size_of_ConcreteBlock[2][1]/2.0;   
  // block #4
  size_of_ConcreteBlock[3][0] = 10.0*ft2mm;
  size_of_ConcreteBlock[3][1] = thick_concrete;
  position_of_ConcreteBlock[3][0] = position_of_ConcreteBlock[2][0] 
                                     - size_of_ConcreteBlock[2][0]/2.0
                                     + size_of_ConcreteBlock[3][0]/2.0;
  position_of_ConcreteBlock[3][1] = position_of_ConcreteBlock[2][1]
                                     - size_of_ConcreteBlock[2][1]/2.0
                                     - size_of_ConcreteBlock[3][1]/2.0;
  // block #5
  size_of_ConcreteBlock[4][0] = 34.0*ft2mm;
  size_of_ConcreteBlock[4][1] = thick_concrete;
  position_of_ConcreteBlock[4][0] = position_of_ConcreteBlock[1][0] 
                                     + size_of_ConcreteBlock[1][0]/2.0
                                     + 6.0*ft2mm
                                     + size_of_ConcreteBlock[4][0]/2.0;                                     
  position_of_ConcreteBlock[4][1] = position_of_ConcreteBlock[1][1];

  // block #6
  size_of_ConcreteBlock[5][0] = thick_concrete;
  size_of_ConcreteBlock[5][1] = 5.0*ft2mm;
  position_of_ConcreteBlock[5][0] = position_of_ConcreteBlock[4][0]
                                     + size_of_ConcreteBlock[4][0]/2.0
                                     - size_of_ConcreteBlock[5][0]/2.0;
  position_of_ConcreteBlock[5][1] = position_of_ConcreteBlock[4][1]
                                     + size_of_ConcreteBlock[4][1]/2.0
                                     + size_of_ConcreteBlock[5][1]/2.0;
  // block #7 
  size_of_ConcreteBlock[6][0] = 52.0*ft2mm;
  size_of_ConcreteBlock[6][1] = thick_concrete;
  position_of_ConcreteBlock[6][0] = position_of_ConcreteBlock[0][0]
                                     - size_of_ConcreteBlock[0][0]/2.0
                                     + size_of_ConcreteBlock[6][0]/2.0; 
  position_of_ConcreteBlock[6][1] = position_of_ConcreteBlock[0][1]
                                     + size_of_ConcreteBlock[0][1]/2.0
                                     + size_of_ConcreteBlock[6][1]/2.0;
  // block #8
  size_of_ConcreteBlock[7][0] = thick_concrete;
  size_of_ConcreteBlock[7][1] = 16.0*ft2mm;
  position_of_ConcreteBlock[7][0] = position_of_ConcreteBlock[6][0]
                                     + size_of_ConcreteBlock[6][0]/2.0
                                     - size_of_ConcreteBlock[7][0]/2.0;
  position_of_ConcreteBlock[7][1] = position_of_ConcreteBlock[6][1]
                                     - size_of_ConcreteBlock[6][1]/2.0
                                     - size_of_ConcreteBlock[7][1]/2.0;
  // block #9
  size_of_ConcreteBlock[8][0] = 12.0*ft2mm;
  size_of_ConcreteBlock[8][1] = thick_concrete;
  position_of_ConcreteBlock[8][0] =  position_of_ConcreteBlock[7][0]
                                     + size_of_ConcreteBlock[7][0]/2.0
                                     - size_of_ConcreteBlock[8][0]/2.0;
  position_of_ConcreteBlock[8][1] =  position_of_ConcreteBlock[7][1]
                                      - size_of_ConcreteBlock[7][1]/2.0
                                      - size_of_ConcreteBlock[8][1]/2.0;

}
//------------------------------------------------
void ELSParameters::detectorplane_parameter_initialize(void)
{
  /* -------------
     Default Value
     ------------- */
  thin_dp  = 1.0;  // unit=mm

  //numDetectorPlane = 198;
 //numDetectorPlane = 58;
  numDetectorPlane = 100;

  //size_of_DetectorPlane[0] = 400.0*1E+3;  //unit = mm
  //size_of_DetectorPlane[1] = 400.0*1E+3;

  //for( int i=0; i<numDetectorPlane; i++ ) {
    //position_of_DetectorPlane[i][0] = 0.0;   // unit = mm
    //position_of_DetectorPlane[i][1] = 0.0;
    //position_of_DetectorPlane[i][2] = 6.0*1E+3 + double(i)*3.0*1E+3;
    //position_of_DetectorPlane[i][2] = 6.0*1E+3 + double(i)*2.0*1E+3;

    //position_of_DetectorPlane[i][0] = beam_injection_position[0];
    //position_of_DetectorPlane[i][1] = beam_injection_position[1];

    /*
    position_of_DetectorPlane[i][2] = 2130.0 + double(i)*100.0;
    if( position_of_DetectorPlane[i][2] <= 4000.0 ){ 
        size_of_DetectorPlane[i][0] = 600.0;  //unit = mm
        size_of_DetectorPlane[i][1] = 600.0;
    }else if( position_of_DetectorPlane[i][2] > 4000.0 && 
              position_of_DetectorPlane[i][2] <= 9000.0 ) {
        size_of_DetectorPlane[i][0] = 4000.0;  //unit = mm
        size_of_DetectorPlane[i][1] = 4000.0;
    }else if( position_of_DetectorPlane[i][2] > 9000.0 ){
        size_of_DetectorPlane[i][0] = 20000.0;  //unit = mm                                                                                
        size_of_DetectorPlane[i][1] = 20000.0;
    }
    size_of_DetectorPlane[i][2] = thin_dp;
    */
  //}

  /*
   Test Detector Plane for check dE/dX by Ionization
   2011.07.20
   */
  //thin_dp  = 5000.0;  // unit=mm, but we should use in 10,20,30,41.1MeV  

  //thin_dp  = 95.0;
  //thin_dp  = 8301.62; 
  //numDetectorPlane = 1;

  //size_of_DetectorPlane[0] = 200.0*1E+3;  //unit = mm
  //size_of_DetectorPlane[1] = 200.0*1E+3;
  //size_of_DetectorPlane[0] = 200.0;  //unit = mm                                                                 
  //size_of_DetectorPlane[1] = 200.0;                          
  //size_of_DetectorPlane[2] = thin_dp;

  //position_of_DetectorPlane[0][0] = 0.0;   // unit = mm
  //position_of_DetectorPlane[0][1] = 0.0;
 
  for( int i=0; i<numDetectorPlane; i++ ) {
       size_of_DetectorPlane[i][0] = 300.0*1E+3;  //unit = mm
       size_of_DetectorPlane[i][1] = 300.0*1E+3;
       size_of_DetectorPlane[i][2] = thin_dp;
  
       position_of_DetectorPlane[i][0] = beam_injection_position[0];
       position_of_DetectorPlane[i][1] = beam_injection_position[1];
       //position_of_DetectorPlane[i][2] = 5001.0 + double(i)*thin_dp + thin_dp/2.0;
       position_of_DetectorPlane[i][2] = 5000.0 + double(i)*2000.0;  // unit=mm
  }

}
//------------------------------------------------
void ELSParameters::inputBeamLineG4TubeList(void)
{
   NumBLComp=0;

    // No.1 EGUN G4Tube[0]     Fragne ICF203-ICF203
    RMinBL[NumBLComp]=RMin_EGUN;
    RMaxBL[NumBLComp]=RMax_EGUN;
    L_BL[NumBLComp]=L_EGUN;
    RMinBL_FF[NumBLComp]=RMin_EGUN_FF;
    RMaxBL_FF[NumBLComp]=RMax_EGUN_FF;
    L_BL_FF[NumBLComp]=L_EGUN_FF;
    RMinBL_BF[NumBLComp]=RMin_EGUN_BF;
    RMaxBL_BF[NumBLComp]=RMax_EGUN_BF;
    L_BL_BF[NumBLComp]=L_EGUN_BF;
    for( int i=0; i<3; i++ ){
     position_BL_FF[NumBLComp][i]=position_EGUN_FF[i];
     position_BL[NumBLComp][i]=position_EGUN[i];
     position_BL_BF[NumBLComp][i]=position_EGUN_BF[i];
     rotation_BL[NumBLComp][i]=rotation_EGUN[i];
    }
    NumBLComp++;

    // No.2 ML         G4Tube     Fragne ICF203-ICF070
    RMinBL[NumBLComp]=RMin_ML;
    RMaxBL[NumBLComp]=RMax_ML;
    L_BL[NumBLComp]=L_ML;
     RMinBL_FF[NumBLComp]=RMin_ML_FF;
     RMaxBL_FF[NumBLComp]=RMax_ML_FF;
     L_BL_FF[NumBLComp]=L_ML_FF;
      RMinBL_BF[NumBLComp]=RMin_ML_BF;
      RMaxBL_BF[NumBLComp]=RMax_ML_BF;
      L_BL_BF[NumBLComp]=L_ML_BF;
    for( int i=0; i<3; i++ ){
      position_BL_FF[NumBLComp][i]=position_ML_FF[i];
      position_BL[NumBLComp][i]=position_ML[i];
      position_BL_BF[NumBLComp][i]=position_ML_BF[i];
      rotation_BL[NumBLComp][i]=rotation_ML[i];
    }
    NumBLComp++;

   // No.3 ICF070-GV  G4Tube        Fragne ICF070-ICF070
    RMinBL[NumBLComp]=RMin_ICF070GV;
    RMaxBL[NumBLComp]=RMax_ICF070GV;
    L_BL[NumBLComp]=L_ICF070GV;
     RMinBL_FF[NumBLComp]=RMin_ICF070GV_FF;
     RMaxBL_FF[NumBLComp]=RMax_ICF070GV_FF;
     L_BL_FF[NumBLComp]=L_ICF070GV_FF;
      RMinBL_BF[NumBLComp]=RMin_ICF070GV_BF;
      RMaxBL_BF[NumBLComp]=RMax_ICF070GV_BF;
      L_BL_BF[NumBLComp]=L_ICF070GV_BF;
    for( int i=0; i<3; i++ ){
      position_BL_FF[NumBLComp][i]=position_ICF070GV_FF[i];
      position_BL[NumBLComp][i]=position_ICF070GV[i];
      position_BL_BF[NumBLComp][i]=position_ICF070GV_BF[i];
      rotation_BL[NumBLComp][i]=rotation_ICF070GV[i];
    }
    NumBLComp++;

   // No.4 ICF070-V85 duct G4Tube   Fragne ICF070-V85
    RMinBL[NumBLComp]=RMin_ICF070V85Duct;
    RMaxBL[NumBLComp]=RMax_ICF070V85Duct;
    L_BL[NumBLComp]=L_ICF070V85Duct;
     RMinBL_FF[NumBLComp]=RMin_ICF070V85Duct_FF;
     RMaxBL_FF[NumBLComp]=RMax_ICF070V85Duct_FF;
     L_BL_FF[NumBLComp]=L_ICF070V85Duct_FF;
      RMinBL_BF[NumBLComp]=RMin_ICF070V85Duct_BF;
      RMaxBL_BF[NumBLComp]=RMax_ICF070V85Duct_BF;
      L_BL_BF[NumBLComp]=L_ICF070V85Duct_BF;
    for( int i=0; i<3; i++ ){
      position_BL_FF[NumBLComp][i]=position_ICF070V85Duct_FF[i];
      position_BL[NumBLComp][i]=position_ICF070V85Duct[i];
      position_BL_BF[NumBLComp][i]=position_ICF070V85Duct_BF[i];
      rotation_BL[NumBLComp][i]=rotation_ICF070V85Duct[i];
    }
    NumBLComp++;
    
    // No.5 CoreMonitor1 G4Tube      Fragne V85-V85
    RMinBL[NumBLComp]=RMin_CM1;
    RMaxBL[NumBLComp]=RMax_CM1;
    L_BL[NumBLComp]=L_CM1;
     RMinBL_FF[NumBLComp]=RMin_CM1_FF;
     RMaxBL_FF[NumBLComp]=RMax_CM1_FF;
     L_BL_FF[NumBLComp]=L_CM1_FF;
      RMinBL_BF[NumBLComp]=RMin_CM1_BF;
      RMaxBL_BF[NumBLComp]=RMax_CM1_BF;
      L_BL_BF[NumBLComp]=L_CM1_BF;
    for( int i=0; i<3; i++ ){
      position_BL_FF[NumBLComp][i]=position_CM1_FF[i];
      position_BL[NumBLComp][i]=position_CM1[i];
      position_BL_BF[NumBLComp][i]=position_CM1_BF[i];
      rotation_BL[NumBLComp][i]=rotation_CM1[i];
    }
    NumBLComp++;

    // No.6 PB+B-Tube G4Tube         Fragne V85-V85
    RMinBL[NumBLComp]=RMin_PBBTube;
    RMaxBL[NumBLComp]=RMax_PBBTube;
    L_BL[NumBLComp]=L_PBBTube;
     RMinBL_FF[NumBLComp]=RMin_PBBTube_FF;
     RMaxBL_FF[NumBLComp]=RMax_PBBTube_FF;
     L_BL_FF[NumBLComp]=L_PBBTube_FF;
      RMinBL_BF[NumBLComp]=RMin_PBBTube_BF;
      RMaxBL_BF[NumBLComp]=RMax_PBBTube_BF;
      L_BL_BF[NumBLComp]=L_PBBTube_BF;
    for( int i=0; i<3; i++ ){
      position_BL_FF[NumBLComp][i]=position_PBBTube_FF[i];
      position_BL[NumBLComp][i]=position_PBBTube[i];
      position_BL_BF[NumBLComp][i]=position_PBBTube_BF[i];
      rotation_BL[NumBLComp][i]=rotation_PBBTube[i];
    }
    NumBLComp++;

    // No.7 CoreMonitor2 G4Tube      Fragne V85-V85
    RMinBL[NumBLComp]=RMin_CM2;
    RMaxBL[NumBLComp]=RMax_CM2;
    L_BL[NumBLComp]=L_CM2;
     RMinBL_FF[NumBLComp]=RMin_CM2_FF;
     RMaxBL_FF[NumBLComp]=RMax_CM2_FF;
     L_BL_FF[NumBLComp]=L_CM2_FF;
      RMinBL_BF[NumBLComp]=RMin_CM2_BF;
      RMaxBL_BF[NumBLComp]=RMax_CM2_BF;
      L_BL_BF[NumBLComp]=L_CM2_BF;
    for( int i=0; i<3; i++ ){
      position_BL_FF[NumBLComp][i]=position_CM2_FF[i];
      position_BL[NumBLComp][i]=position_CM2[i];
      position_BL_BF[NumBLComp][i]=position_CM2_BF[i];
      rotation_BL[NumBLComp][i]=rotation_CM2[i];
    }
    NumBLComp++;
 
    // No.8 2m-Tube G4Tube           Fragne V85-V85
    RMinBL[NumBLComp]=RMin_2mTube;
    RMaxBL[NumBLComp]=RMax_2mTube;
    L_BL[NumBLComp]=L_2mTube;
     RMinBL_FF[NumBLComp]=RMin_2mTube_FF;
     RMaxBL_FF[NumBLComp]=RMax_2mTube_FF;
     L_BL_FF[NumBLComp]=L_2mTube_FF;
      RMinBL_BF[NumBLComp]=RMin_2mTube_BF;
      RMaxBL_BF[NumBLComp]=RMax_2mTube_BF;
      L_BL_BF[NumBLComp]=L_2mTube_BF;
    for( int i=0; i<3; i++ ){
      position_BL_FF[NumBLComp][i]=position_2mTube_FF[i];
      position_BL[NumBLComp][i]=position_2mTube[i];
      position_BL_BF[NumBLComp][i]=position_2mTube_BF[i];
      rotation_BL[NumBLComp][i]=rotation_2mTube[i];
    }
    NumBLComp++;

    // No.9 Drift Duct G4Tube        Fragne V85-V85
    RMinBL[NumBLComp]=RMin_V85V85Duct;
    RMaxBL[NumBLComp]=RMax_V85V85Duct;
    L_BL[NumBLComp]=L_V85V85Duct;
     RMinBL_FF[NumBLComp]=RMin_V85V85Duct_FF;
     RMaxBL_FF[NumBLComp]=RMax_V85V85Duct_FF;
     L_BL_FF[NumBLComp]=L_V85V85Duct_FF;
      RMinBL_BF[NumBLComp]=RMin_V85V85Duct_BF;
      RMaxBL_BF[NumBLComp]=RMax_V85V85Duct_BF;
      L_BL_BF[NumBLComp]=L_V85V85Duct_BF;
    for( int i=0; i<3; i++ ){
      position_BL_FF[NumBLComp][i]=position_V85V85Duct_FF[i];
      position_BL[NumBLComp][i]=position_V85V85Duct[i];
      position_BL_BF[NumBLComp][i]=position_V85V85Duct_BF[i];
      rotation_BL[NumBLComp][i]=rotation_V85V85Duct[i];
    }
    NumBLComp++;

    // No.10 V85-GV G4Tube           Fragne V85-V85
    RMinBL[NumBLComp]=RMin_V85GV;
    RMaxBL[NumBLComp]=RMax_V85GV;
    L_BL[NumBLComp]=L_V85GV;
     RMinBL_FF[NumBLComp]=RMin_V85GV_FF;
     RMaxBL_FF[NumBLComp]=RMax_V85GV_FF;
     L_BL_FF[NumBLComp]=L_V85GV_FF;
      RMinBL_BF[NumBLComp]=RMin_V85GV_BF;
      RMaxBL_BF[NumBLComp]=RMax_V85GV_BF;
      L_BL_BF[NumBLComp]=L_V85GV_BF;
    for( int i=0; i<3; i++ ){
      position_BL_FF[NumBLComp][i]=position_V85GV_FF[i];
      position_BL[NumBLComp][i]=position_V85GV[i];
      position_BL_BF[NumBLComp][i]=position_V85GV_BF[i];
      rotation_BL[NumBLComp][i]=rotation_V85GV[i];
    }
    NumBLComp++;

    // No.11 QM G4Tube               Fragne V85-V85
    RMinBL[NumBLComp]=RMin_QM;
    RMaxBL[NumBLComp]=RMax_QM;
    L_BL[NumBLComp]=L_QM;
     RMinBL_FF[NumBLComp]=RMin_QM_FF;
     RMaxBL_FF[NumBLComp]=RMax_QM_FF;
     L_BL_FF[NumBLComp]=L_QM_FF;
      RMinBL_BF[NumBLComp]=RMin_QM_BF;
      RMaxBL_BF[NumBLComp]=RMax_QM_BF;
      L_BL_BF[NumBLComp]=L_QM_BF;
    for( int i=0; i<3; i++ ){
      position_BL_FF[NumBLComp][i]=position_QM_FF[i];
      position_BL[NumBLComp][i]=position_QM[i];
      position_BL_BF[NumBLComp][i]=position_QM_BF[i];
      rotation_BL[NumBLComp][i]=rotation_QM[i];
    }
    NumBLComp++;

    // No.13 CoreMonitor3 G4Tube     Fragne V85-V85
    RMinBL[NumBLComp]=RMin_CM3;
    RMaxBL[NumBLComp]=RMax_CM3;
    L_BL[NumBLComp]=L_CM3;
     RMinBL_FF[NumBLComp]=RMin_CM3_FF;
     RMaxBL_FF[NumBLComp]=RMax_CM3_FF;
     L_BL_FF[NumBLComp]=L_CM3_FF;
      RMinBL_BF[NumBLComp]=RMin_CM3_BF;
      RMaxBL_BF[NumBLComp]=RMax_CM3_BF;
      L_BL_BF[NumBLComp]=L_CM3_BF;
    for( int i=0; i<3; i++ ){
      position_BL_FF[NumBLComp][i]=position_CM3_FF[i];
      position_BL[NumBLComp][i]=position_CM3[i];
      position_BL_BF[NumBLComp][i]=position_CM3_BF[i];
      rotation_BL[NumBLComp][i]=rotation_CM3[i];
    }
    NumBLComp++;

    // No.15 V85-GV2 G4Tube           Fragne V85-V85
    RMinBL[NumBLComp]=RMin_V85GV2;
    RMaxBL[NumBLComp]=RMax_V85GV2;
    L_BL[NumBLComp]=L_V85GV2;
     RMinBL_FF[NumBLComp]=RMin_V85GV2_FF;
     RMaxBL_FF[NumBLComp]=RMax_V85GV2_FF;
     L_BL_FF[NumBLComp]=L_V85GV2_FF;
      RMinBL_BF[NumBLComp]=RMin_V85GV2_BF;
      RMaxBL_BF[NumBLComp]=RMax_V85GV2_BF;
      L_BL_BF[NumBLComp]=L_V85GV2_BF;
    for( int i=0; i<3; i++ ){
      position_BL_FF[NumBLComp][i]=position_V85GV2_FF[i];
      position_BL[NumBLComp][i]=position_V85GV2[i];
      position_BL_BF[NumBLComp][i]=position_V85GV2_BF[i];
      rotation_BL[NumBLComp][i]=rotation_V85GV2[i];
    }
    NumBLComp++;

    // No.16 T-Duct G4Tube           Fragne V85-V85
    RMinBL[NumBLComp]=RMin_TDuct;
    RMaxBL[NumBLComp]=RMax_TDuct;
    L_BL[NumBLComp]=L_TDuct;
     RMinBL_FF[NumBLComp]=RMin_TDuct_FF;
     RMaxBL_FF[NumBLComp]=RMax_TDuct_FF;
     L_BL_FF[NumBLComp]=L_TDuct_FF;
      RMinBL_BF[NumBLComp]=RMin_TDuct_BF;
      RMaxBL_BF[NumBLComp]=RMax_TDuct_BF;
      L_BL_BF[NumBLComp]=L_TDuct_BF;
    for( int i=0; i<3; i++ ){
      position_BL_FF[NumBLComp][i]=position_TDuct_FF[i];
      position_BL[NumBLComp][i]=position_TDuct[i];
      position_BL_BF[NumBLComp][i]=position_TDuct_BF[i];
      rotation_BL[NumBLComp][i]=rotation_TDuct[i];
    }
    NumBLComp++;


    // No.18 DriftTube G4Tube        Fragne V85-V85
    RMinBL[NumBLComp]=RMin_V85V85Duct2;
    RMaxBL[NumBLComp]=RMax_V85V85Duct2;
    L_BL[NumBLComp]=L_V85V85Duct2;
     RMinBL_FF[NumBLComp]=RMin_V85V85Duct2_FF;
     RMaxBL_FF[NumBLComp]=RMax_V85V85Duct2_FF;
     L_BL_FF[NumBLComp]=L_V85V85Duct2_FF;
      RMinBL_BF[NumBLComp]=RMin_V85V85Duct2_BF;
      RMaxBL_BF[NumBLComp]=RMax_V85V85Duct2_BF;
      L_BL_BF[NumBLComp]=L_V85V85Duct2_BF;
    for( int i=0; i<3; i++ ){
      position_BL_FF[NumBLComp][i]=position_V85V85Duct2_FF[i];
      position_BL[NumBLComp][i]=position_V85V85Duct2[i];
      position_BL_BF[NumBLComp][i]=position_V85V85Duct2_BF[i];
      rotation_BL[NumBLComp][i]=rotation_V85V85Duct2[i];
    }
    NumBLComp++;

    // No.20 Drift Duct G4Tube       Fragne V85-V85
    RMinBL[NumBLComp]=RMin_V85V85Duct3;
    RMaxBL[NumBLComp]=RMax_V85V85Duct3;
    L_BL[NumBLComp]=L_V85V85Duct3;
     RMinBL_FF[NumBLComp]=RMin_V85V85Duct3_FF;
     RMaxBL_FF[NumBLComp]=RMax_V85V85Duct3_FF;
     L_BL_FF[NumBLComp]=L_V85V85Duct3_FF;
      RMinBL_BF[NumBLComp]=RMin_V85V85Duct3_BF;
      RMaxBL_BF[NumBLComp]=RMax_V85V85Duct3_BF;
      L_BL_BF[NumBLComp]=L_V85V85Duct3_BF;
    for( int i=0; i<3; i++ ){
      position_BL_FF[NumBLComp][i]=position_V85V85Duct3_FF[i];
      position_BL[NumBLComp][i]=position_V85V85Duct3[i];
      position_BL_BF[NumBLComp][i]=position_V85V85Duct3_BF[i];
      rotation_BL[NumBLComp][i]=rotation_V85V85Duct3[i];
    }
    NumBLComp++;

    // No.22 CoreMonitor4 G4Tube     Fragne V85-V85
    RMinBL[NumBLComp]=RMin_CM4;
    RMaxBL[NumBLComp]=RMax_CM4;
    L_BL[NumBLComp]=L_CM4;
     RMinBL_FF[NumBLComp]=RMin_CM4_FF;
     RMaxBL_FF[NumBLComp]=RMax_CM4_FF;
     L_BL_FF[NumBLComp]=L_CM4_FF;
      RMinBL_BF[NumBLComp]=RMin_CM4_BF;
      RMaxBL_BF[NumBLComp]=RMax_CM4_BF;
      L_BL_BF[NumBLComp]=L_CM4_BF;
    for( int i=0; i<3; i++ ){
      position_BL_FF[NumBLComp][i]=position_CM4_FF[i];
      position_BL[NumBLComp][i]=position_CM4[i];
      position_BL_BF[NumBLComp][i]=position_CM4_BF[i];
      rotation_BL[NumBLComp][i]=rotation_CM4[i]; 
    }
    NumBLComp++;

}
//====================================================
void ELSParameters::dump_ELS_geometory(void)
{
  cout << "Concrete Pad" << endl;
  cout << " Size " 
       << " " << size_of_concretepad[0] 
       << " " << size_of_concretepad[1]
       << " " << size_of_concretepad[2] 
       << endl;
  cout << " Position "
       << " " << position_of_concretepad[0]
       << " " << position_of_concretepad[1]
       << " " << position_of_concretepad[2]
       << endl;
   
  cout << "ELS Container (Outer)" << endl;
  cout << " Size "
       << " " << outer_size_of_ELScontainer[0]
       << " " << outer_size_of_ELScontainer[1]
       << " " << outer_size_of_ELScontainer[2]
       << endl;
  cout << " Position "
       << " " << position_of_ELS_outer[0]
       << " " << position_of_ELS_outer[1]
       << " " << position_of_ELS_outer[2]
       << endl;

  cout << "ELS Container (Inner)" << endl;
  cout << " Size "
       << " " << inner_size_of_ELScontainer[0]
       << " " << inner_size_of_ELScontainer[1]
       << " " << inner_size_of_ELScontainer[2]
       << endl;
  cout << " Position "
       << " " << position_of_ELS_inner[0]
       << " " << position_of_ELS_inner[1]
       << " " << position_of_ELS_inner[2]
       << endl;

  cout << "ELS Cover Box" << endl;
  cout << " Position "
       << " " << position_of_coverbox_outer[0]
       << " " << position_of_coverbox_outer[1]
       << " " << position_of_coverbox_outer[2] 
       << endl;

  for( int i=0; i<3; i++ ){
    cout << "Base of " << i << endl;
    cout << " Size "
	 << " " << size_of_base[i][0]
	 << " " << size_of_base[i][1]
	 << " " << size_of_base[i][2]
	 << endl;
    cout << " Position "
	 << " " << position_of_base[i][0]
	 << " " << position_of_base[i][1]
	 << " " << position_of_base[i][2]
	 << endl;
  }

  for( int i=0; i<NumBLComp; i++ ){
    cout << "Beam Line " << i << endl;
  
    cout << " Position "
         << " " <<  position_BL_FF[i][0] 
         << " " <<  position_BL[i][0] 
         << " " <<  position_BL_BF[i][0] 
	 << endl;
  }

  
  


}


