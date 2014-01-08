/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*%												   %*/
/*% Functions: ComputeHorizon and ComputeViewfactors	 			                   %*/
/*%												   %*/
/*% N. Helbig                                                                                      %*/
/*% Last changes: 2008 									           %*/
/*%												   %*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "debug.h"

#ifdef _PAROC_
#include "EnergyBalance.ph"
#else
#include "EnergyBalance.h"
#define rprintf printf
#endif

#include "Constants.h"

int EnergyBalance::ComputeHorizon( int ee, 
				   int ff, 
				   int gg, 
				   int hh )
{

 /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
 /*%                                                                                                %*/
 /*% Last changes: 2008                                                                             %*/
 /*%                                                                                                %*/
 /*% "ray-tracing"-method for the horizon calculation:                                              %*/
 /*% for every grid cell the horizon is determined by comparing the horizon angles of every grid    %*/
 /*% cell along a beam to every border grid cell                                                    %*/
 /*%                                                                                                %*/
 /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
                                 
 double t1, t0;                 // clock variables for measuring the time needed in the main ComputeHorizon 
 double delay;                  // for the delay of calling the clock function
 int a, b;                      // control variables
 double number;                 // control variable
 long int count_vf = 0;         // counts number of View factors(i,j,a,b) > 0

 // for the horizon detection:
 double azi;                    // azimuth angle of the viewing grid cell beam to the border grid cells
 double flatdist;               // flat distance for calculating the azimuth angle to the border cell
 double tMaxX, tMaxY;           // maximum ray length components that are necessary to reach the next grid cell
 double tDeltaX, tDeltaY;       // maximum ray length components to cross the width of one grid cell under azi 
 int stepsx, stepsy;            // ariables that indicate if i,j have to be in- or decremented to reach target 
 double bx, by, bz;             // x, y, z - distances between the surfaces
 bool visible;                  // for finding the visibility/horizon angle
 int ho_a, ho_b;                // x,y of the horizon per viewed grid point
 double view_angle;             // elevation angle for target grid cell along the beam
 double testview_angle;         // angle for terrain height along the beam of a grid cell
 
 // for the view factor calculation:
 double d_distij;               // distance of the respective plane to the origin
 double iza, izb, izc, izd;     // z-values of the corner vectors around grid center i+0.5,j+0.5
 double dist;                   // actual 3D - center distance between surfaces
 double area_ij;                // large unsplitted surface areas

 bool sensor_reached = false; 
 double backup_z;               // backup for the real i,j -> z
 double backup_sx;              // backup for the real i,j -> sx
 double backup_sy;              // backup for the real i,j -> sy
 
 t0 = clock();
 t1 = clock();
 delay = t1 - t0;
 
 // for loops over all i,j - grid cells 
 for ( int i = 0; i < dimx; i++ )
 {
     for ( int j = 0; j < dimy; j++ )
     {
         // at model start from 'InitializeViewFactors'
	 if ( model_start )
	 {
	    vf_t[i][j] = 0.;	
	 }
	 // if "!store_vf" or "calc_suncontrol" or !model_start (for the sensor) from 'ComputeRadiationBalance'
	 else if ( !model_start || !store_vf || calc_suncontrol )
	 {
	    i = ee;
            j = ff;
	 }
         else
	 {
	 }	 

       	 // first the specific sensor plane is computed
	 if ( i == sensor.i && j == sensor.j && !sensor_reached )
	 {
            sensor_reached = true;		 
            // backup for real i,j properties
	    backup_z  = z[i][j];
	    backup_sx = sx[i][j];
	    backup_sy = sy[i][j];
	    // sensor properties are temporarily given to the real i,j
	    z[i][j]  = sensor.z;
            sx[i][j] = sensor.sx;
            sy[i][j] = sensor.sy;
	 }
         // afterwards the real i,j plane is computed	 
	 else if ( sensor_reached )
	 {
	    sensor_reached = false;
	    i = sensor.i;
	    j = sensor.j;
	    
	    // since this grid cell coordinates are computed twice vf_t has to be set to zero again
	    if ( model_start )
	    {
	       vf_t[i][j] = 0.;
	    }
	    else
	    {
	    }
	    // backup is given back to the "real" patch
	    z[i][j]  = backup_z;
	    sx[i][j] = backup_sx;
	    sy[i][j] = backup_sy;
	 }
	 else
	 {
	 }		 
	
         // for loop over border grid cells
	 for ( int m = 0; m < dimx; m++ )
       	 {
	     for ( int t = 0; t < dimy; t++ )
	     {
	         // horizon is only calculated by looking to every border grid cell BUT
		 // a mutual visibility test is needed between selected or all grid cells (depending on store_vf) as well as
		 // for the sun-or-shadow indicator a visibility test to a border grid cell (sun position)
		 if ( !model_start || !store_vf || calc_suncontrol )
		 {
		    m = gg;
	            t = hh;
		    vf[0][0] = 0;
		 }	    
                 
                 if ( m != i || m != i || t != j || t != j )
	         {
	
	            if ( calc_suncontrol )
		    {
		       // special border treatment of ij to assure that for border cells the horizon is, by looking along the
		       // border, the farthest that grid cell can view
		       // this is done by giving that cell the appropriate edge cell to determine the horizon for that
			                                                      
                       if ( (j == 0 && t == 0) || (j == dimy - 1 && t == dimy - 1) )
		       {
			  if ( m < i )
			  {    
			     switch ( j )
	                     {
	                         case 0:
	                                m = 0;
				        t = 0;
				 break;
				 default:
				        m = 0;
				        t = dimy - 1;
				 break;
		 	     }
			  }    
			  else
			  {
			     switch ( j )
			     {
			         case 0:
			                m = dimx - 1;
		                        t = 0;
		                 break; 
	                         default:
	                                m = dimx - 1;
	                                t = dimy - 1;				       
				 break;
			     }  
			  }  
		       } // end of if condition : upper and lower border 
			  
                       if ( (i == 0 && m == 0) || (i == dimx - 1 && m == dimx -1) )		
		       {
			  if ( t < j )
			  {
			     switch ( i )
			     {
			         case 0:	    
				        m = 0;
				        t = 0;
				 break;	   
				 default:
			                 m = dimx - 1;
			                 t = 0;
			         break;
		             }
	                  }  
                          else
			  {
	                     switch ( i )
	                     {
	                         case 0:
	                                m = 0;
                                        t = dimy - 1;
                                 break;
                                 default:
                                        m = dimx - 1;
                                        t = dimy - 1;
                                 break;
			     }
			  }  				
		       } // end of if condition : left and right border	 

		    } // end of if special border treatment for sun-or-shadow detection

	            ho_a = m;
	            ho_b = t;   

		    // the view angle is derived from the target grid cell:
		    // than a more efficient request for detecting if the grid cell is 
		    // in sun or not resp. visible is possible by sampling if a view angle
		    // on the way to the target grid cell is larger in that way in most cases 
		    // one does not to go all way to the border of the model domain
		    
		    // distance of the ij-plane to the origin
		    // grid cell z-value is shifted to center
		    d_distij = sx[i][j] * (i + 0.5) * cellsize + sy[i][j] * (j + 0.5) * cellsize
		              + z[i][j] * 1.;

		    // z-values of the edge vectors at the forced to be rectangular area ij:
		    // (i,j), (i+1,j), (i,j+1) and (i+1,j+1)
		    // they are forced on the ij-plane, z-component of the normal vector is 1 (cf. 'CalculateAziSlope')
                    iza = (d_distij - sx[i][j] * (i       * cellsize)
		                    - sy[i][j] * (j       * cellsize)) / 1.;
                    izb = (d_distij - sx[i][j] * ((i + 1) * cellsize)
		                    - sy[i][j] * (j       * cellsize)) / 1.;
	            izc = (d_distij - sx[i][j] * (i       * cellsize)
		                    - sy[i][j] * ((j + 1) * cellsize)) / 1.;
		    izd = (d_distij - sx[i][j] * ((i + 1) * cellsize)
		                    - sy[i][j] * ((j + 1) * cellsize)) / 1.;
                     
		    // bx = dx (m), by = dy (m), bz = dz (m)
                    // to get the actual viewing beam vector
                    // the following is valid also for inclined surfaces
                    // by forcing the corner z-values on the plane with sx, sy, sz
		    bx = (m - i) * cellsize;
		    by = (t - j) * cellsize;
		    bz = z[m][t] - z[i][j]; 
	
		    dist = sqrt( (bx * bx) + (by * by) + (bz * bz) );

		    // view angle is in case of sun-or-shadow detection the solar elevation angle
		    if ( calc_suncontrol )
		    {
		       view_angle = sin( solar.elev );
		    }
		    // and in the other case the elevation angle of the target grid cell
		    else
		    {
		       view_angle = bz / dist;
		    }

		    // round the view angle to avoid that internal rounding or inaccuracies
		    // prevent that the lower testview angle and view angle comparison go wrong
		    view_angle *= pow(10, 4);
		    if ( view_angle >= 0 )
		    {
	               number = floor( view_angle + 0.5 );
		    }
		    else
		    {
		       number = ceil( view_angle - 0.5 );
		    }
		    view_angle = number / pow(10, 4);
		                                      
                    // calculation of pseudo azimuth angle for every border grid cell:
		    // from south counterclockwise
		    flatdist = sqrt( (i - m) * (i - m) + (j - t) * (j - t) );
		    azi = acos( -(t - j) / flatdist );

		    if ( m - i > 0 )
		    {
		    }
		    else
		    {
		       azi = 2. * PI - azi;
		    }
                    
		    // variables that indicate if i,j have to be in- or decremented to cross border cells
		    if ( m - i > 0 )
		    {
		       stepsx = 1;	    
		    }
		    else if ( m - i < 0 )
		    {
		       stepsx = -1;	    
		    }
		    else
		    {
		       stepsx = 0;	    
		    }
		    
		    if ( t - j > 0 )
		    {
		       stepsy = 1;
		    }
		    else if ( t - j < 0 )
		    {
		       stepsy = -1;
		    }
		    else
		    {
		       stepsy = 0;    
		    }

		    if ( fabs( stepsx ) > 0 )
		    {
		       // the ray length components that are necessary to reach the next grid cell border
		       // in x = i direction starting from the center of ij:
		       tMaxX =   (0.5 * cellsize) / sin( azi );
		       // the maximum ray length component x to cross the width of one grid cell (voxel) 
		       // under a certain azimuth angle
		       tDeltaX =   cellsize / sin( azi );
		    }
		    else
		    {
		       tMaxX = 0;	    
		       tDeltaX = 0;
		    }
		    
		    if ( fabs( stepsy ) > 0 ) 
		    {
		       // the ray length components that are necessary to reach the next grid cell border
		       // in y = j direction starting from the center of ij:	    
		       tMaxY = - (0.5 * cellsize) / cos( azi );
		       // the maximum ray length component y to cross the width of one grid cell (voxel)
		       // under a certain azimuth angle
		       tDeltaY = - cellsize / cos( azi );
		    }
		    else
		    {
		       tMaxY = 0;
		       tDeltaY = 0;
		    }
	
	            // round tMaxX/tMaxY/tDeltaX/tDeltaY to avoid that internal rounding or inaccuracies
		    // prevent that the comparison between them go wrong
		         	    
                    tMaxX *= pow(10, 4);
		    if ( tMaxX >= 0 )
                    {
	               number = floor( tMaxX + 0.5 );
	            }
	            else
	            {
                       number = ceil( tMaxX - 0.5 );
                    }
                    tMaxX = number / pow(10, 4);

		    tDeltaX *= pow(10, 4);
                    if ( tDeltaX >= 0 )
                    {
                       number = floor( tDeltaX + 0.5 );
                    }
                    else
                    {
                       number = ceil( tDeltaX - 0.5 );
                    }
                    tDeltaX = number / pow(10, 4);
										
		    tMaxY *= pow(10, 4);
                    if ( tMaxY >= 0 )
                    {
                       number = floor( tMaxY + 0.5 );
                    }
		    else
		    {
		       number = ceil( tMaxY - 0.5 );
		    }
		    tMaxY = number / pow(10, 4);

                    tDeltaY *= pow(10, 4);
                    if ( tDeltaY >= 0 )
                    {
                       number = floor( tDeltaY + 0.5 );
                    }
                    else
                    {
                       number = ceil( tDeltaY - 0.5 );
                    }
                    tDeltaY = number / pow(10, 4);
		    
		    a = i;
		    b = j; 
		    
		    visible = true;
		    
                    // do while loop over view angle comparison until target is reached
		    do
	            {
		      // visibility test according to Amanatides(1987) for ray tracing	   
		       
	              if ( (fabs( tMaxX ) < fabs( tMaxY ) && fabs( stepsx ) > 0) || tMaxY == 0 )
		      {
		         a = a + stepsx;	
		         tMaxX = tMaxX + tDeltaX;
		      }
		      else if ( (fabs( tMaxY ) < fabs( tMaxX ) && fabs( stepsy ) > 0) || tMaxX == 0 )
		      {
			 b = b + stepsy;
			 tMaxY = tMaxY + tDeltaY;
		      }
		      else // for tMaxX==tMaxY (45, 135, 225, 315 degree azimuth)
		      {
			 a = a + stepsx;
		         tMaxX = tMaxX + tDeltaX;
		         b = b + stepsy;
		         tMaxY = tMaxY + tDeltaY;	 
		      }

		      if ( a < 0 || a > dimx - 1 || b < 0 || b > dimy - 1 )
		      {  
		         // beam vector components lie outside of the model region
			 break;                        
		      }
                   
		      // bx = dx (m), by = dy (m), bz = dz (m) 
		      // to get the actual viewing beam vector
		      // the following is valid also for inclined surfaces
		      // by forcing the corner z-values on the plane with sx, sy, sz 
		      bx = (a - i) * cellsize;
		      by = (b - j) * cellsize;
                      bz = z[a][b] - z[i][j];  
		     
		      dist = sqrt( (bx * bx) + (by * by) + (bz * bz) );
		      
		      testview_angle = bz / dist;
		      // round the testing view angle to avoid that internal rounding or inaccuracies
		      // prevent that the following testing view angle and view angle comparison go wrong
		      testview_angle *= pow(10, 4);
		      if ( testview_angle >= 0 )
		      {
		         number = floor( testview_angle + 0.5 );
		      }
		      else
		      {
		         number =  ceil( testview_angle - 0.5 );
		      }
		      testview_angle = number / pow(10, 4);
											   
		      if ( dist > 0. )
		      {
		         // note that since the plane can be inclined after substructuring some small elements might 
		         // see each other although bz was zero, i.e. comparing only the center coordinates
		         // is not sufficient
			 if ( testview_angle > view_angle )
		         {
		            ho_a = a;
		            ho_b = b;		  
			    visible = false;       // terrain higher inbetween
		         }
		         else
		         {
		         }
		      }	     
		      // the compared grid cells are identical (dist = 0) 
		      else    
		      {  
		      }
		      
		      // a,b is the target border grid cell m,t
		      if ( a == m && b == t )
		      {
		         break;
		      }
		    
	            } while ( visible ); // end of do while loop over view angle comparison until target grid cell is reached
			    
                    // higher terrain was found inbetween
	            if ( !visible )
		    {
		    }  
		    else
		    {
		       /* THE VIEW FACTOR CALCULATION: */
			  
                       // all view factors are stored or are not stored for 
		       if ( !calc_suncontrol )
		       {
                          ComputeViewfactors ( i, j, a, b, area_ij, count_vf, sensor_reached );
		       }
		       // only for the horizon of the sun-or-shadow detection
		       else
		       {
		       }
		       
                    } // end of if two grid cells or sun projected cell are mutually visible
		
	            // the actual calculated horizon coordinates are stored but only needed for sun-or-shadow detection	       
		    ho[0][0] = ho_a;     
		    ho[0][1] = ho_b;
		       
		    if ( sensor_reached )
		    {
		       sensor.hox = ho_a;
		       sensor.hoy = ho_b;
		    }
	            else
		    {
		    }
		       
		    if ( !model_start || !store_vf || calc_suncontrol ) 
		    {
		       if ( sensor_reached )
                       {
		          m = dimx;
			  t = dimy;
		       }
		       else
		       {
                          m = dimx;
                          t = dimy;
		          i = dimx;
		          j = dimy;
		       }
		     }
		     else
		     {
		     }
                    
                  } // end of if condition : only for m,t that are not equal i,j        
	          else
	          {
	             // the grid cell coordinates of ij as horizon to itself are stored as horizon coordinates 
		     // but only needed if sun-or-shadow is detected
		     ho[0][0] = i;
		     ho[0][1] = j;

		     if ( sensor_reached )
		     {
		        sensor.hox = i;
	                sensor.hoy = j;
		     }	       
		     else
		     {
		     }
		    
		     if ( !model_start || !store_vf || calc_suncontrol )
		     {	       
		        if ( sensor_reached )
		        {
		 	   m = dimx;
			   t = dimy; 
		        }
		        else
		        {
			   i = dimx;
			   j = dimy;
			   m = dimx;
			   t = dimy;
		        }
		     }
                  } // end of else condition : only for m,t that are equal i,j
	      } // end of for over t  
          } // end of for over m (all grid cells) 

		   
	  if ( model_start ) 
	  {
             if ( vf_t[i][j] > 1.0 ) 
	     {
	        rprintf( "\n at i %u j %u sum terrain view factors = %lf > 1;", i, j, vf_t[i][j] );
	        fflush( stdout );	
	     }
	    
             // sky view factor calculation from the total sum of terrain view factors at grid cell
	     sky_vf[i][j] = 1.0 - vf_t[i][j];

	     if ( sensor_reached )
	     {
	        sensor.sky_vf = sky_vf[i][j];
	     }
	    
	     // total terrain view factor multiplied by the actual area for a different usage in the radiation exchanging part
       	     vf_t[i][j] *= area_ij;
          }
          else if ( !store_vf && !calc_suncontrol )
          {
             if ( vf[0][0] > 0.9 )
             {
                rprintf( "\n at i %u j %u to a %u b %u terrain view factor = %lf > 1;", ee, ff, gg, hh, vf[0][0] );
                rprintf( "\n better check if a lower area substructuring threshold for" );
                rprintf( "\n the view factor computation is necessary" );
                rprintf( "\n :in your ebalance/EBParam.ini" );
                fflush( stdout );
             }
          }
          else
          {
          }
      }  // end of for j loop over all j grid cells  
  } // end of for i loop over all i grid cells

  if ( model_start )
  {
      t1 = clock();
      rprintf("\nend of viewfactor computation: time needed %2.10lf sum of Fij>0 %4.2lf\% \n",
  		                   (t1 - t0 - delay) / CLOCKS_PER_SEC,
 	     ((double)count_vf * 100.)/(dimx * dimy * dimx * dimy));
  }
  else
  {
  }
  return(EXIT_SUCCESS);

} // end of ComputeHorizon

int EnergyBalance::ComputeViewfactors( int i, 
		                       int j, 
				       int a, 
				       int b, 
				       double &area_ij,
	                               long int &count_vf,
                                       bool sensor_reached )
{
	
 /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
 /*%                                                                                                    %*/
 /*% Last changes: 2008                                                                                 %*/ 
 /*%                                                                                                    %*/
 /*% calculation of the View Factor that describes the fraction of the total power that an              %*/
 /*% area A_i emits (diffuse) and that arrives on A_j (dimensionless geometric factor)                  %*/
 /*% 													%*/
 /*% for close grid cells:                                                                              %*/
 /*% area to area view factor calculation for every pair of grid cells that see each other using the    %*/
 /*% substructuring method                                                                              %*/
 /*% for distant grid cells:                                                                            %*/
 /*% point to area view factor calculation for every pair of grid cells that see each other from the    %*/
 /*% center of the respective area                                                                      %*/
 /*%                                                                                                    %*/
 /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/ 	
 
 int ij, ab;                    // control variables;	
 double bx, by, bz;             // distances between the viewers large and the viewed large surface
 double d_distij, d_distab;     // distance of the respective plane to the origin
 double iza, izb, izc, izd;     // z-values of the vectors in the center of the 3 areas around grid point ij
 double aza, azb, azc, azd;     // z-values of the vectors in the center of the 3 areas around grid point ab
 double lengthi, lengthj;       // actual grid cell dimensions of ij considering slope angle
 double lengtha, lengthb;       // actual grid cell dimensions of ab considering slope angle
 int spliti, splitj;            // splitting factors for the respective grid cell dimension
 int splita, splitb;            // splitting factors for the respective grid cell dimension
 double area_ab;                // large unsplitted surface areas
 double area_uv, area_kl;       // small surface areas because of splitting uv in ij-plane and kl in ab-plane
 double x_uv, y_uv, z_uv;       // vector components of the small cells in the ij-plane
 double x_kl, y_kl, z_kl;       // vector components of the small cells in the ab-plane
 double dist;                   // actual centercenter  distance between every small pair of areas
 double cosangle_ij;            // angle between normal to [i,j] and beam from [i,j]
 double cosangle_ab;            // angle between normal to [a,b] and beam from [i,j]
 double add_vf;                 // for summing up the small view factors in the substructured grid cell

 ij = j * dimx + i;
 ab = b * dimx + a;

 // grid cell z-value is shifted to center coordinate
 // distance of the ij-plane to the origin
 d_distij = sx[i][j] * (i + 0.5) * cellsize + sy[i][j] * (j + 0.5) * cellsize
           + z[i][j]  * 1.;

 // distance of the ab-plane to the origin
 d_distab = sx[a][b] * (a + 0.5) * cellsize + sy[a][b] * (b + 0.5) * cellsize
           + z[a][b]  * 1.;

 // z-values of the edge vectors at the forced to be rectangular area ij: (i,j), (i+1,j), (i,j+1) and (i+1,j+1)
 // they are forced on the ij-plane around a centered coordinate, z-component of the normal vector is 1 ('CalculateAziSlope')
 iza = (d_distij - sx[i][j] * (i       * cellsize)
		 - sy[i][j] * (j       * cellsize)) / 1.;
 izb = (d_distij - sx[i][j] * ((i + 1) * cellsize)
                 - sy[i][j] * (j       * cellsize)) / 1.;
 izc = (d_distij - sx[i][j] * (i       * cellsize)
                 - sy[i][j] * ((j + 1) * cellsize)) / 1.;
 izd = (d_distij - sx[i][j] * ((i + 1) * cellsize)
		 - sy[i][j] * ((j + 1) * cellsize)) / 1.;

 // and likewise for the z-values of the edge vectors at area ab
 aza = (d_distab - sx[a][b] * (a       * cellsize)
                 - sy[a][b] * (b       * cellsize)) / 1.;
 azb = (d_distab - sx[a][b] * ((a + 1) * cellsize)
                 - sy[a][b] * (b       * cellsize)) / 1.;
 azc = (d_distab - sx[a][b] * (a       * cellsize)
		 - sy[a][b] * ((b + 1) * cellsize)) / 1.;
 azd = (d_distab - sx[a][b] * ((a + 1) * cellsize)
		 - sy[a][b] * ((b + 1) * cellsize)) / 1.;
 
 // find the longest side of area ij
 lengthi = sqrt( ((i + 1) * cellsize - i * cellsize) *
                 ((i + 1) * cellsize - i * cellsize)
               + (j       * cellsize - j * cellsize) *
                 (j       * cellsize - j * cellsize)
               + (izb - iza) * (izb - iza) );
 lengthj = sqrt( (i       * cellsize - i * cellsize) *
                 (i       * cellsize - i * cellsize)
               + ((j + 1) * cellsize - j * cellsize) *
                 ((j + 1) * cellsize - j * cellsize)
               + (izc - iza) * (izc - iza) );
 
 // find the longest side of area ab 
 lengtha = sqrt( ((a + 1) * cellsize - a * cellsize) *
                 ((a + 1) * cellsize - a * cellsize)
               + (b       * cellsize - b * cellsize) *
                 (b       * cellsize - b * cellsize)
               + (azb - aza) * (azb - aza) );
 lengthb = sqrt( (a       * cellsize - a * cellsize) *
                 (a       * cellsize - a * cellsize)
               + ((b + 1) * cellsize - b * cellsize) *
                 ((b + 1) * cellsize - b * cellsize)
               + (azc - aza) * (azc - aza) );
 
 // surface area = |n|, assumed that the 3D-area's are parallelograms:
 area_ij = sqrt( sx[i][j] * sx[i][j] + sy[i][j] * sy[i][j] + 1. )
                 * cellsize * cellsize;
 area_ab = sqrt( sx[a][b] * sx[a][b] + sy[a][b] * sy[a][b] + 1. )
                 * cellsize * cellsize;

 // the following is valid also for inclined surfaces as the planes are enlarged 
 // by forcing the corner z-values on the plane with sx, sy, sz at the lower left
 // corner of the grid cell
 bx = (a - i) * cellsize;               // bx = dx (m) going from [a,b] to [i,j]
 by = (b - j) * cellsize;               // by = dy (m) going from [a,b] to [i,j]
 bz = z[a][b] - z[i][j];                // bz = dz (m) going from [a,b] to [i,j]
 
 // center distance between the large surfaces
 dist = sqrt( (bx * bx) + (by * by) + (bz * bz) );
 
 // substructuring if the side length of area ij divided by the distance exceeds 
 // the preassumed value 'sub_crit'
 // McCluney (1994) assumes that the maximum side length should be lower than 10% 
 // of the distance as a rough criteria until when a finite source can still be 
 // treated as a point source (no substructuring is required)
 if ( lengthi / dist >= sub_crit ||
      lengthj / dist >= sub_crit ||
      lengtha / dist >= sub_crit ||
      lengthb / dist >= sub_crit )
 {
    // computing the necessary splitting factors for all border lengths
    if ( lengthi / dist >= sub_crit )
    {
       spliti = (int)((lengthi / (dist * sub_crit)) + 0.5);
    }                                	    
    else
    {	    
       spliti = 1;
    }
    
    if ( lengthj / dist >= sub_crit )
    {
       splitj = (int)((lengthj / (dist * sub_crit)) + 0.5);
    }
    else
    {
       splitj = 1;
    }
    
    if ( spliti < 1 )
    { 
        spliti = 1;
    }
    
    if ( splitj < 1 )
    {
        splitj = 1;
    }
    
    // find the largest splitting factor for both sides
    if ( spliti >= splitj )
    {							                                   
       splitj = spliti;
    }
    else
    {
       spliti = splitj;
    }
    
    // likewise for the visible area ab
    if ( lengtha / dist >= sub_crit )  
    {
       splita = (int)((lengtha / (dist * sub_crit)) + 0.5);	    
    }
    else
    {
       splita = 1;
    }
    
    if ( lengthb / dist >= sub_crit )    
    {
       splitb = (int)((lengthb / (dist * sub_crit)) + 0.5);
    }
    else
    {
       splitb = 1;
    }

    if ( splita < 1 )
    {
       splita = 1;
    }
    
    if ( splitb < 1 )
    {
       splitb = 1;
    }
    
    if ( splita >= splitb )
    {
       splitb = splita;
    }
    else
    {
       splita = splitb;
    }

    // small surface areas in the subdivided areas ij and ab
    area_uv = area_ij / (spliti * splitj);
    area_kl = area_ab / (splita * splitb);
    
    add_vf = 0.;                  // every small view factor needs to be calculated and summed   
    for ( int u = 0; u < spliti; u++ )
    {
        for ( int v = 0; v < splitj; v++ )
        {
            // find the coordinate values of the small cell uv in the plane ij
	    // vec(a) + vec(b)-vec(a) * (u + 0.5) / spliti
	    //        + vec(c)-vec(a) * (v + 0.5) / splitj
	    x_uv = i   +  ((i + 1) - i) * (u + 0.5) * (1. / spliti)
                              + (i - i) * (v + 0.5) * (1. / splitj);
            y_uv = j   +        (j - j) * (u + 0.5) * (1. / spliti)
	               +  ((j + 1) - j) * (v + 0.5) * (1. / splitj);
	    z_uv = iza + (izb - iza)    * (u + 0.5) * (1. / spliti)
	               + (izc - iza)    * (v + 0.5) * (1. / splitj);
	    
	    for ( int k = 0; k < splita; k++ )
	    {
                for ( int l = 0; l < splitb; l++ )
		{
                    // find the coordinate values of the small cell kl in the plane ab
		    // vec(a) + vec(b)-vec(a) * (k + 0.5) / splita
		    //        + vec(c)-vec(a) * (l + 0.5) / splitb
                    x_kl = a   + ((a + 1) - a) * (k + 0.5) * (1. / splita)
                                     + (a - a) * (l + 0.5) * (1. / splitb);
                    y_kl = b   +       (b - b) * (k + 0.5) * (1. / splita)
		               + ((b + 1) - b) * (l + 0.5) * (1. / splitb);
		    z_kl = aza + (azb - aza)   * (k + 0.5) * (1. / splita)
		               + (azc - aza)   * (l + 0.5) * (1. / splitb);
									
		    // new distance between small cell area_uv and small area_kl
		    dist = sqrt( (x_kl - x_uv) * (x_kl - x_uv) * cellsize * cellsize
                               + (y_kl - y_uv) * (y_kl - y_uv) * cellsize * cellsize
                               + (z_kl - z_uv) * (z_kl - z_uv) );
		    
		    // cosangle_ij: n * (b) / (|n| * |b|)
		    // the angle between normal to area_ij and beam from area_uv
		    cosangle_ij = (sx[i][j] * (x_kl - x_uv) * cellsize
		   	         + sy[i][j] * (y_kl - y_uv) * cellsize + 1. * (z_kl - z_uv))
			         / (sqrt( sx[i][j] * sx[i][j] +
			                  sy[i][j] * sy[i][j] + 1. ) * dist);
		    
		    // cosangle_ab: n * (-b) / (|n| * |b|)
		    // the angle between normal to area_ab and beam from area_uv
		    cosangle_ab = (sx[a][b] * (x_uv - x_kl) * cellsize
                                 + sy[a][b] * (y_uv - y_kl) * cellsize + 1. * (z_uv - z_kl))
                                 / (sqrt( sx[a][b] * sx[a][b] +
                                          sy[a][b] * sy[a][b] + 1. ) * dist);
		    
		    if ( cosangle_ij > 0. && cosangle_ab > 0. )
		    {
		        add_vf += (cosangle_ij * cosangle_ab * area_uv * area_kl) / (PI * dist * dist);
		    }
		}    
            }   		    
	}    	    
    }   	    
    if ( model_start && store_vf )
    {
       // dividing through the whole surface area according to the view factor formula
       vf[ij][ab] = add_vf / area_ij;
       if ( add_vf > 0. ) count_vf++;
       vf_t[i][j] += vf[ij][ab];
    }
    else if ( model_start && !store_vf )
    {
       vf[0][0] =  add_vf / area_ij;
       if ( add_vf > 0. ) count_vf++; 
       vf_t[i][j] += vf[0][0];
    } 
    else
    {
       vf[0][0] =  add_vf / area_ij;	   
       if ( sensor_reached ) 
       {
	  sensor.vf = add_vf / area_ij;
       }
    }
 } // end of if condition : if a length of an area > 10 % of the distance
 else
 {
    add_vf = 0;	 
    
    // cosangle_ij: n * (b) / (|n| * |b|)
    // the angle between normal to area_ij and beam from area_ij
    cosangle_ij = (sx[i][j] * bx + sy[i][j] * by + 1. * bz) /
                  (sqrt( (sx[i][j] * sx[i][j]) + (sy[i][j] * sy[i][j]) + 1. )
                 * sqrt( dist * dist ));
     
    // cosangle_ab: n * (-b) / (|n| * |b|)
    // the angle between normal to area_ab and beam from area_ij
    cosangle_ab = (sx[a][b] * (-bx) + sy[a][b] * (-by) + 1. * (-bz)) /
                  (sqrt( (sx[a][b] * sx[a][b]) + (sy[a][b] * sy[a][b]) + 1. )
                 * sqrt( dist * dist ));
     
    if ( cosangle_ij > 0. && cosangle_ab > 0. )
    {  
       add_vf = (cosangle_ij * cosangle_ab * area_ab) /
                (PI * dist * dist);	       
       
       if ( model_start && store_vf )
       {	       
          vf[ij][ab] = add_vf;
	  if ( vf[ij][ab] > 0. ) count_vf++;
	  vf_t[i][j] += vf[ij][ab];
       }
       else if ( model_start && !store_vf )
       {
	  vf[0][0] =  add_vf;
	  if ( vf[0][0] > 0. ) count_vf++;
	  vf_t[i][j] += vf[0][0];
       }
       else
       {
	  vf[0][0] =  add_vf;
	  if ( sensor_reached )
	  {
	     sensor.vf = add_vf;
	  }	  
       }
    }
 }
 return EXIT_SUCCESS;
} // end of ComputeViewfactors
