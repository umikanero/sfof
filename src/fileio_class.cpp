/*Class for file input and output*/

#include "fileio_class.hpp"

void Fileio::split (const std::string &str, std::vector<std::string> &tokens, 
		    const std::string &delimiter) {
  //! Function to split a string read in from a file into columns.
  std::string::size_type lastPos, pos;      
  lastPos = str.find_first_not_of(delimiter, 0); /* pos = find first "non-delimiter" */
  pos = str.find_first_of(delimiter, lastPos);
  while (std::string::npos != pos || std::string::npos != lastPos){
    tokens.push_back(str.substr(lastPos, pos-lastPos)); /* Found a token, add it to the vector. */     
    lastPos = str.find_first_not_of(delimiter, pos); /* Skip delimiters. Note the "not_of" */         
    pos=str.find_first_of(delimiter, lastPos); /* Find next "non-delimiter" */            
  }
}

void Fileio::read_ascii (const std::string &fname, const std::string &mode, 
			 std::vector<Galaxy> &gals) { 
  //! Function to read in an ASCII file and store the contents in a vector of Galaxy instances.
  int count = 0; /* line count */
  std::string line; /* line string */
  std::vector<std::string> cols; /* column vector */
  std::ifstream read_file(fname.c_str()); /* open file */
  while(!read_file.eof()) { /* while not the end of the file */
    std::getline(read_file, line); /* read each line */
    if(line.length() >= 1 && 
       line.find("#") == std::string::npos) { /* skip empty lines and lines starting with # */
      split(line, cols, " "); /* split line into columns */
      if(mode == "spec") {
	Galaxy spec_gal(count, cols[0], atof(cols[1].c_str()), atof(cols[2].c_str()), 
			atof(cols[3].c_str())); /* intialise spec Galaxy */
	gals.push_back(spec_gal); /* store spec Galaxy instance in vector */
      }
      else {
	Galaxy phot_gal(count, cols[0], atof(cols[1].c_str()), atof(cols[2].c_str()), 
			atof(cols[3].c_str()), atof(cols[4].c_str())); /* intialise phot Galaxy */
	gals.push_back(phot_gal); /* store phot Galaxy instance in vector */
      }
      cols.clear(); /* clear column vector */
      count++; /* increase line count */
    }
  }
  read_file.close(); /* close file */
}

void Fileio::read_fits (const std::string &fname, const std::string &mode,
			std::vector<Galaxy> &gals) { 
  //! Function to read in an FITS file and store the contents in a vector of Galaxy instances.
  fitsfile *fptr;
  int n_cols, status=0, hdunum, hdutype, anynul, typecode;
  long n_rows, repeat, width;
  char *val, nullstr[] = "*";
  std::vector<std::string> cols; /* column vector */
  val = new char[1000];
  if(!fits_open_file(&fptr, fname.c_str(), READONLY, &status)){ /* open FITS file */
    if(fits_get_hdu_num(fptr, &hdunum) == 1) fits_movabs_hdu(fptr, 2, &hdutype, &status);
    else fits_get_hdu_type(fptr, &hdutype, &status);
    if(hdutype == IMAGE_HDU) /* check if the FITS file contains an image */
      std::cout<<"Error: this program only displays tables, not images.\n"<<std::flush;
    else{
      fits_get_num_rows(fptr, &n_rows, &status); /* get number of rows in FITS table */
      fits_get_num_cols(fptr, &n_cols, &status); /* get number of columns in FITS table */
      for(int i = 1; i <= n_rows && !status; i++){
	for(int j = 1; j <= n_cols; j++){
	  fits_get_coltype(fptr, j, &typecode, &repeat, &width, &status);
	  for(int k = 1; k <= repeat; k++){
	    if(fits_read_col_str(fptr, j, i, k, 1, nullstr, &val, &anynul, &status)) break;
	    cols.push_back(val);
	  }
	}
	if(mode == "spec") {
	  Galaxy spec_gal(i - 1, cols[0], atof(cols[1].c_str()), atof(cols[2].c_str()), 
			  atof(cols[3].c_str())); /* intialise spec Galaxy */
	  gals.push_back(spec_gal); /* store spec Galaxy instance in vector */
	}
	else {
	  Galaxy phot_gal(i - 1, cols[0], atof(cols[1].c_str()), atof(cols[2].c_str()), 
			  atof(cols[3].c_str()), atof(cols[4].c_str())); /* intialise phot Galaxy */
	  gals.push_back(phot_gal); /* store phot Galaxy instance in vector */
	}
	cols.clear(); /* clear column vector */
      }      
    }
  }
}

void Fileio::output_file_names (const std::string &fname, const std::string &mode, 
				const std::string &output, double link_r, double link_z) {
  //! Function to set up output file names.
  std::stringstream cluster_file_stream, member_file_stream;
  if(output == "ascii") {
    cluster_file_stream<<fname<<"_clusters_"<<link_r<<"_"<<link_z<<"_"<<mode<<".dat";
    member_file_stream<<fname<<"_members_"<<link_r<<"_"<<link_z<<"_"<<mode<<".dat"; 
  }
  else if(output == "fits") {
    cluster_file_stream<<fname<<"_clusters_"<<link_r<<"_"<<link_z<<"_"<<mode<<".fits";
    member_file_stream<<fname<<"_members_"<<link_r<<"_"<<link_z<<"_"<<mode<<".fits"; 
  }
  cluster_file_name = cluster_file_stream.str();
  member_file_name = member_file_stream.str();
}

void Fileio::write_ascii (const std::vector<Cluster> &cluster_list) {
  //! Funtion to output Cluster intances to an ASCII file.
  std::ofstream write_clusters (cluster_file_name.c_str());
  std::ofstream write_members (member_file_name.c_str());
  for(int i = 0; i < cluster_list.size(); i++) {
    write_clusters<<std::fixed<<std::setprecision(3);
    write_clusters<<std::setw(6)<<cluster_list[i].num<<" ";
    write_clusters<<std::setw(6)<<cluster_list[i].ngal<<" ";
    write_clusters<<std::setw(7)<<cluster_list[i].ra<<" ";
    write_clusters<<std::setw(7)<<std::showpos<<cluster_list[i].dec<<" ";
    write_clusters<<std::setw(5)<<std::noshowpos<<cluster_list[i].z<<" ";
    write_clusters<<std::setw(5)<<cluster_list[i].size<<" ";
    write_clusters<<std::setw(7)<<cluster_list[i].area<<"\n";
    for(int j = 0; j < cluster_list[i].ngal; j++) {
      write_members<<std::fixed<<std::setprecision(3);
      write_members<<std::setw(6)<<cluster_list[i].num<<" ";
      write_members<<std::setw(6)<<cluster_list[i].ngal<<" ";
      write_members<<std::setw(12)<<cluster_list[i].mem[j].id<<" ";
      write_members<<std::setw(7)<<cluster_list[i].mem[j].ra<<" ";
      write_members<<std::setw(7)<<std::showpos<<cluster_list[i].mem[j].dec<<" ";
      write_members<<std::setw(5)<<std::noshowpos<<cluster_list[i].mem[j].z<<"\n";
    }
  }
  write_clusters.close();
  write_members.close();
}

void Fileio::write_fits (const std::vector<Cluster> &cluster_list) {
  //! Funtion to output Cluster intances to an FITS file.
  fitsfile *fptr1, *fptr2;
  int tint, current_pos = 0, status = 0;
  double tdouble;
  const int cluster_fields = 7, member_fields = 6;
  char *cluster_types[cluster_fields], *member_types[member_fields];
  char *cluster_forms[cluster_fields], *member_forms[member_fields], *tstring;
  cluster_types[0] = const_cast<char *>("NUM");   cluster_forms[0] = const_cast<char *>("J");
  cluster_types[1] = const_cast<char *>("NGAL");  cluster_forms[1] = const_cast<char *>("J");
  cluster_types[2] = const_cast<char *>("RA");    cluster_forms[2] = const_cast<char *>("E");
  cluster_types[3] = const_cast<char *>("DEC");   cluster_forms[3] = const_cast<char *>("E");
  cluster_types[4] = const_cast<char *>("Z");     cluster_forms[4] = const_cast<char *>("E");
  cluster_types[5] = const_cast<char *>("SIZE");  cluster_forms[5] = const_cast<char *>("E");
  cluster_types[6] = const_cast<char *>("AREA");  cluster_forms[6] = const_cast<char *>("E");
  member_types[0] = const_cast<char *>("C_NUM");  member_forms[0] = const_cast<char *>("J");
  member_types[1] = const_cast<char *>("C_NGAL"); member_forms[1] = const_cast<char *>("J");
  member_types[2] = const_cast<char *>("G_ID");   member_forms[2] = const_cast<char *>("12A");
  member_types[3] = const_cast<char *>("G_RA");   member_forms[3] = const_cast<char *>("E");
  member_types[4] = const_cast<char *>("G_DEC");  member_forms[4] = const_cast<char *>("E");
  member_types[5] = const_cast<char *>("G_Z");    member_forms[5] = const_cast<char *>("E");
  fits_create_file(&fptr1, cluster_file_name.c_str(), &status); /*create new FITS file*/
  fits_create_file(&fptr2, member_file_name.c_str(), &status); /*create new FITS file*/
  if(status != 0){
    std::cout<<"Error! Cannot create FITS file.\n"<<std::flush;
    exit(-1);
  }
  fits_create_tbl(fptr1, BINARY_TBL, 0, cluster_fields, cluster_types, cluster_forms, NULL, NULL, &status); /*create new FITS table*/ 
  fits_create_tbl(fptr2, BINARY_TBL, 0, member_fields, member_types, member_forms, NULL, NULL, &status); /*create new FITS table*/ 
  for(int i = 1; i <= cluster_list.size(); i++) {
    tint = cluster_list[i - 1].num;
    fits_write_col(fptr1, TUINT, 1, i, 1, 1, &tint, &status);
    tint = cluster_list[i - 1].ngal;
    fits_write_col(fptr1, TUINT, 2, i, 1, 1, &tint, &status);
    tdouble = cluster_list[i - 1].ra;
    fits_write_col(fptr1, TDOUBLE, 3, i, 1, 1, &tdouble, &status);
    tdouble = cluster_list[i - 1].dec;
    fits_write_col(fptr1, TDOUBLE, 4, i, 1, 1, &tdouble, &status);
    tdouble = cluster_list[i - 1].z;
    fits_write_col(fptr1, TDOUBLE, 5, i, 1, 1, &tdouble, &status);
    tdouble = cluster_list[i - 1].size;
    fits_write_col(fptr1, TDOUBLE, 6, i, 1, 1, &tdouble, &status);
    tdouble = cluster_list[i - 1].area;
    fits_write_col(fptr1, TDOUBLE, 7, i, 1, 1, &tdouble, &status);
    for(int j = 1; j <= cluster_list[i - 1].ngal; j++) {
      current_pos++;
      tint = cluster_list[i - 1].num;
      fits_write_col(fptr2, TUINT, 1, current_pos, 1, 1, &tint, &status);
      tint = cluster_list[i - 1].ngal;
      fits_write_col(fptr2, TUINT, 2, current_pos, 1, 1, &tint, &status);
      tstring = const_cast<char *>(cluster_list[i - 1].mem[j - 1].id.c_str());
      fits_write_col(fptr2, TSTRING, 3, current_pos, 1, 1, &tstring, &status); 
      tdouble = cluster_list[i - 1].mem[j - 1].ra;
      fits_write_col(fptr2, TDOUBLE, 4, current_pos, 1, 1, &tdouble, &status);
      tdouble = cluster_list[i - 1].mem[j - 1].dec;
      fits_write_col(fptr2, TDOUBLE, 5, current_pos, 1, 1, &tdouble, &status);
      tdouble = cluster_list[i - 1].mem[j - 1].z;
      fits_write_col(fptr2, TDOUBLE, 6, current_pos, 1, 1, &tdouble, &status);
    }
  }
  fits_close_file(fptr1, &status); /*close FITS file*/
  fits_close_file(fptr2, &status); /*close FITS file*/
}