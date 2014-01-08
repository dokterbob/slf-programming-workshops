double RadiationField::getHorizon(const size_t& ix1, const size_t& iy1, const double& bearing) const
{
	bool horizon_found = false;
	size_t nb_cells = 1;
	const double cell_alt = dem.grid2D(ix1, iy1);
	double horizon_tan_angle = 0.;
	const double alpha = bearing*to_rad;

	if(ix1==0 || ix1==dem_dimx-1 || iy1==0 || iy1==dimy-1) return 0.; //a border cell is not shadded

	while(!horizon_found) {
		nb_cells++;
		const int ix2 = (int)ix1 + (int)round( ((double)nb_cells)*sin(alpha) ); //alpha is a bearing
		const int iy2 = (int)iy1 + (int)round( ((double)nb_cells)*cos(alpha) ); //alpha is a bearing

		if(ix2<=0 || ix2>=(signed)dem_dimx-1 || iy2<=0 || iy2>=(signed)dimy-1) break; //we are out of the dem

		const double new_altitude = dem.grid2D(ix2, iy2);
		if(new_altitude==mio::IOUtils::nodata) break; //we stop at nodata cells

		const double DeltaH = new_altitude - cell_alt;
		const double distance = sqrt( (double)( (ix2-ix1)*(ix2-ix1) + (iy2-iy1)*(iy2-iy1)) ) * cellsize;
		const double tan_angle = DeltaH/distance;
		if(tan_angle>horizon_tan_angle) horizon_tan_angle = tan_angle;

		if(distance>max_shade_distance) horizon_found=true; //maximum lookup distance reached
	}

	return horizon_tan_angle;
}