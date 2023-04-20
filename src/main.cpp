#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
//#include <rsf.h>
#include <assert.h>
#include "headers.h"
#include <gptl.h>


int main(int argc, char** argv){
    double t_start, t_end;
    double t_start0;

    int size, dims[3];

    int reorder = 1; // true
    int periods[3] = {0, 0, 0}; // false

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


	PARAMS params;
	GRID grid = { 0 };
		
	MPI_NEIGHBOR mpiNeigbor = { 0 };

	getParams( &params );
    // Get the parameters from the command line
    get_params(argc, argv);
    

    NT = (int)(TMAX / DT);
    RNT = (int)((NT + 8 - 1) / 8 ) * 8;

    // Create Seismic Wave Communication domain
    dims[0] = PX; dims[1] = PY; dims[2] = PZ;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, reorder, &SWMPI_COMM);
    MPI_Cart_coords(SWMPI_COMM, this_rank, 3, thisid);
    MPI_Cart_shift(SWMPI_COMM, 0, 1, &neigxid[0], &neigxid[1]);
    MPI_Cart_shift(SWMPI_COMM, 1, 1, &neigyid[0], &neigyid[1]);
    MPI_Cart_shift(SWMPI_COMM, 2, 1, &neigzid[0], &neigzid[1]);
    
    MPI_COORD thisMPICoord = { thisid[0], thisid[1], thisid[2] };
    MPI_Comm comm_cart = SWMPI_COMM;

    int color_xy = this_rank%PZ;
    MPI_Comm sub_comm_xy;
    MPI_Comm_split(MPI_COMM_WORLD, color_xy, this_rank, &sub_comm_xy);
    int sub_rank_xy;
    MPI_Comm_rank(sub_comm_xy, &sub_rank_xy);
    int sub_size_xy;
    MPI_Comm_size(sub_comm_xy, &sub_size_xy);

    int color_x = sub_rank_xy/Gather_SIZE_Y;
    MPI_Comm sub_comm_x;
    MPI_Comm_split(sub_comm_xy, color_x, sub_rank_xy, &sub_comm_x);
    MPI_Comm_rank(sub_comm_x, &sub_rank_x);
    int sub_size_x;
    MPI_Comm_size(sub_comm_x, &sub_size_x);

    int color_y_all = sub_rank_xy%PY;
    MPI_Comm sub_comm_y_all;
    MPI_Comm_split(sub_comm_xy, color_y_all, sub_rank_xy, &sub_comm_y_all);
    int sub_rank_y_all;
    MPI_Comm_rank(sub_comm_y_all, &sub_rank_y_all);
    int sub_size_y_all;
    MPI_Comm_size(sub_comm_y_all, &sub_size_y_all);

    int color_y = sub_rank_y_all/Gather_SIZE_X;
    MPI_Comm sub_comm_y;
    MPI_Comm_split(sub_comm_y_all, color_y, sub_rank_y_all, &sub_comm_y);
    MPI_Comm_rank(sub_comm_y, &sub_rank_y);
    int sub_size_y;
    MPI_Comm_size(sub_comm_y, &sub_size_y);



    int ret;
    ret = GPTLinitialize ();               // Initialize GPTL 
    ret = Gstart("main");

	init_grid( params, &grid, thisMPICoord );
    ni = grid.nx; nx = ni + 6; // include ghost points
    nj = grid.ny; ny = nj + 6;
    nk = grid.nz; nz = nk + 6;
    ni1 = 3; ni2 = ni + 3;
    nj1 = 3; nj2 = nj + 3;
    nk1 = 3; nk2 = nk + 3;
    nx1 = 0; nx2 = nx;
    ny1 = 0; ny2 = ny;
    nz1 = 0; nz2 = nz;
    ngi1 = grid.frontNX + ni1; ngi2 = grid.frontNX + ni2;
    ngj1 = grid.frontNY + nj1; ngj2 = grid.frontNY + nj2;
    ngk1 = grid.frontNZ + nk1; ngk2 = grid.frontNZ + nk2;

    //for pml
    isx1 = 0, isx2 = 0, isy1 = 0, isy2 = 0, isz1 = 0, isz2 = 0;

    //for 128B aligned padding
    lnz = (nz * sizeof(float) + 127) / 128 * 128 / sizeof(float);
    lnz = nz;
    if(this_rank == 0) printf("Aligned padding lnz = %d\n", lnz);

    isMaster(this_rank);
    if(this_rank == 0) masternode = 1;
#ifdef FreeSurface
    if(neigzid[1] == MPI_PROC_NULL) {isz2 = 0; freenode = 1;}
#else
    if(neigzid[1] == MPI_PROC_NULL) {absnode = 1; isz2 = 1;}
#endif

    if(neigxid[0] == MPI_PROC_NULL) {absnode = 1; isx1 = 1;}
    if(neigxid[1] == MPI_PROC_NULL) {absnode = 1; isx2 = 1;}
    if(neigyid[0] == MPI_PROC_NULL) {absnode = 1; isy1 = 1;}
    if(neigyid[1] == MPI_PROC_NULL) {absnode = 1; isy2 = 1;}
    if(neigzid[0] == MPI_PROC_NULL) {absnode = 1; isz1 = 1;}


    // set station info
    int stationNum;
	STATION station;
	stationNum = readStationIndex( grid );
	if (stationNum > 0) {
		allocStation( &station, stationNum, NT);
		initStationIndex( grid, station); 
	}
    
    PGV pgv;
    if(freenode == 1) allocatePGV(grid, &pgv);

    // get receiver parameters
    get_nrec();
    if(masternode) printf("debug: nrec = %d, RNT = %d\n", nrec, RNT);
    int *rec_index = (int*)malloc(sizeof(int) * 4 * nrec);
    float *R = (float*)malloc(sizeof(float)*RNT*nrec*3);
    int i;
    for(i = 0; i < 4*nrec; i++) rec_index[i] = 0;
    for(i = 0; i < RNT*nrec*3; i++) R[i] = 0;
    // Load_Receiver(rec_index, this_rank);
    // if(masternode) for(i=0; i<nrec; i++) 
    //     printf("rec_index: this_rank = %d, r_index = %d, ix = %d, iy = %d\n",
    //             rec_index[i*4+0], rec_index[i*4+1], rec_index[i*4+2], rec_index[i*4+3]);

    // Print some parameters to the screen
    if(masternode) paras_msg(); 

    swmpi_datatype();
    t_start0 = MPI_Wtime();

    // Alloc and construct C
    float *C = Alloc_coord();



    if(masternode) printf("constructing coord ...\n");
    //construct_coord(C);

	constructCoord( grid, thisMPICoord, C );

    MPI_Barrier(MPI_COMM_WORLD);
    if(masternode) printf("Preprocess Terrain ...\n");
	preprocessTerrain( params, comm_cart, thisMPICoord, grid, C  );	
    MPI_Barrier(MPI_COMM_WORLD);
    if(masternode) printf("# start to read source ...\n");
    
    // set data to write out
    SLICE_DATA sliceData, sliceDataCpu;
	SLICE slice = { 0 };
	locateSlice(params, grid, &slice);
	allocSliceData(grid, slice, &sliceData);
   
    long long * srcIndex, * srcSlicePos;
    long long pointSrc;
    
    SOURCE_FILE_INPUT src_in = { 0 };
    MOMENT_RATE momentRate;
    dealSource( params, grid, C, thisMPICoord, &srcIndex, &srcSlicePos, &momentRate, &src_in);
    MPI_Barrier(MPI_COMM_WORLD);
    if (masternode) printf( "alloc moment rate slice\n" );
    MOMENT_RATE momentRateSlice;
	if ( 0 != src_in.npts  )
	{
        // verifySrcIndex( thisMPICoord, srcIndex, src_in.npts );
	     momentRateSlice = allocMomentRateSlice(src_in);
	}

    MPI_Barrier(MPI_COMM_WORLD);

    // Alloc Qs damp and wave
    float *Qs = Alloc_attenu();
    float *M = Alloc_metric();
    //float *D = Alloc_media();
    // float *hW = Alloc_wave();
    // float *mW = Alloc_wave();
    // float *tW = Alloc_wave();
    // float *W = Alloc_wave();

    float *hW_8 = Alloc_wave_8();
    float *mW_8 = Alloc_wave_8();
    float *tW_8 = Alloc_wave_8();
    float *W_8 = Alloc_wave_8();
    float *hW_1 = Alloc_wave_1();
    float *mW_1 = Alloc_wave_1();
    float *tW_1 = Alloc_wave_1();
    float *W_1 = Alloc_wave_1();

    float *pml;
    struct aux Aux;
#ifdef usePML 
    pml = Alloc_pml();
    Alloc_aux(&Aux, isx1, isx2, isy1, isy2, isz1, isz2, this_rank);
    abs_init(pml);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    Grid3D damp = NULL;
    damp = Alloc3D(nx, ny, lnz);
    wp_yzs0 = Alloc1D(ny * lnz * WSIZE * LenFDL);
    wp_yzr0 = Alloc1D(ny * lnz * WSIZE * LenFDL);
    wp_yzs1 = Alloc1D(ny * lnz * WSIZE * LenFDL);
    wp_yzr1 = Alloc1D(ny * lnz * WSIZE * LenFDL);
    wp_xzs0 = Alloc1D(nx * lnz * WSIZE * LenFDL);
    wp_xzr0 = Alloc1D(nx * lnz * WSIZE * LenFDL);
    wp_xzs1 = Alloc1D(nx * lnz * WSIZE * LenFDL);
    wp_xzr1 = Alloc1D(nx * lnz * WSIZE * LenFDL);
    wp_xys0 = Alloc1D(nx * ny * WSIZE * LenFDL);
    wp_xyr0 = Alloc1D(nx * ny * WSIZE * LenFDL);
    wp_xys1 = Alloc1D(nx * ny * WSIZE * LenFDL);
    wp_xyr1 = Alloc1D(nx * ny * WSIZE * LenFDL);
    float *matVx2Vz = Alloc1D(nx*ny*9);
    float *matVy2Vz = Alloc1D(nx*ny*9);
    MPI_Barrier(MPI_COMM_WORLD);
    if(masternode) printf("Alloc finished\n");
    if(masternode) printf("extend coord ...\n");
    extend_coord(C, CSIZE); 
    if(masternode) printf("exchange coord ...\n");
    MPI_Barrier(MPI_COMM_WORLD);
    coord_exchange(C);
    if(masternode) printf("finished exchange coord ...\n");
    MPI_Barrier(MPI_COMM_WORLD);

    float range[2];
    cal_range_steph(C, range);
    float hmin = range[0];
    float hmax = range[1];
    float hmin_global, hmax_global;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&hmin, &hmin_global, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&hmax, &hmax_global, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);


    if(masternode) printf("calculate metric ...\n");
    cal_metric(C,M);

    if(masternode) printf("extend metric ...\n");
    extend_Symm_array(M, MSIZE);

    if(masternode) printf("exchange metric ...\n");
    MPI_Barrier(MPI_COMM_WORLD);
    metric_exchange(M);
    MPI_Barrier(MPI_COMM_WORLD);
    if(masternode) printf("finished exchange metric ...\n");

    if(masternode) printf("constructing media ...\n");
    //modified by cbw
    //if(icoord) construct_media(C, M);
    constructMedium(grid, M);

    //construct_media(C,M);
    dealWeisenShenMedium( params, grid, M, C, thisMPICoord );
    //if ( masternode )
    //    printf( "==================M = %f\n", M[MSIZE * 10000 + 11] );
    MPI_Barrier(MPI_COMM_WORLD);
    // data2D_Model_out( thisMPICoord, params, grid, C, M, slice, sliceData, SAMP, SAMPLE_SIZE, sub_rank_x, sub_comm_x, sub_rank_y, sub_comm_y, Gather );
    MPI_Barrier(MPI_COMM_WORLD);

    if(masternode) printf("finished constructing media ...\n");


    float media_range[6];
    cal_range_media(M, media_range);
    if(masternode) printf("finished cal_range_media ...\n");
    float  vp_min = media_range[0];
    float  vp_max = media_range[1];
    float  vs_min = media_range[2];
    float  vs_max = media_range[3];
    float rho_min = media_range[4];
    float rho_max = media_range[5];
    float vp_min_global, vs_min_global, rho_min_global;
    float vp_max_global, vs_max_global, rho_max_global;
    MPI_Barrier(MPI_COMM_WORLD);
    //if(masternode) printf("=========================split line=========================\n");
    MPI_Reduce(&vp_min, &vp_min_global, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&vp_max, &vp_max_global, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&vs_min, &vs_min_global, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&vs_max, &vs_max_global, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&rho_min, &rho_min_global, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&rho_max, &rho_max_global, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
    //printf("range h = %10.2e ~ %10.2e, in this_rank %d %d %d\n", hmin, hmax, thisid[0], thisid[1], thi    sid[2]);

    float dtmax = 1.3f * hmin_global / vp_max_global;

    if(masternode){
        printf("global range of h   = %10.2e ~ %10.2e m\n", hmin_global, hmax_global);
        printf("global range of vp  = %10.2e ~ %10.2e m/s\n", vp_min_global, vp_max_global);
        printf("global range of vs  = %10.2e ~ %10.2e m/s\n", vs_min_global, vs_max_global);
        printf("global range of rho = %10.2e ~ %10.2e kg/m^3\n", rho_min_global, rho_max_global);
        if(DT < dtmax){
            printf("DT = %10.2e < dtmax = %10.2e (sec) satisfy stability condition, OK\n", DT, dtmax);
        }else{
            printf("Serious Error: DT = %10.2e > dtmax = %10.2e (sec) do not satisfy stability condition, ABORT!\n", DT, dtmax);
            MPI_Abort(MPI_COMM_WORLD, 110);
        }
    }

    /*
     * *************************************************************
     * Read Source Info
     * *************************************************************
     */

    coef_surface(M, matVx2Vz, matVy2Vz, pml);
    MPI_Barrier(MPI_COMM_WORLD);

    if(masternode) printf("calculate cerjan coefs ...\n");
    abs_exp(damp); // Caculate cerjan coefs

    if(masternode)
        printf("time loop beginning ...\n");

    int it = 0; 
    //float t = 0.0;
    char snapname[400];

    char filenamebasex[100];
    char filenamebasey[100];
    char filenamebasez[100];
    char filenamebaseeta[100];
    char filenamebaseep[100];
    sprintf(filenamebasex, "%s/SX", OUT);
    sprintf(filenamebasey, "%s/SY", OUT);
    sprintf(filenamebasez, "%s/SZ", OUT);
    sprintf(filenamebaseeta, "%s/Eta", OUT);
    sprintf(filenamebaseep, "%s/EP", OUT);

    int rec_NX, rec_NY, rec_NZ;
    int rec_nbgx, rec_nbgy, rec_nbgz;
    int rec_nedx, rec_nedy, rec_nedz;
    int rec_nxt, rec_nyt, rec_nzt;
    MPI_Offset displacement;

    // same for each processor:
    if (NEDX == -1) NEDX = NX;
    if (NEDY == -1) NEDY = NY;
    if (NEDZ == -1) NEDZ = NZ;
    // make NED's a record point
    // for instance if NBGX:NSKPX:NEDX = 1:3:9
    // then we have 1,4,7 but NEDX=7 is better
    NEDX = NEDX - (NEDX - NBGX) % NSKPX;
    NEDY = NEDY - (NEDY - NBGY) % NSKPY;
    NEDZ = NEDZ - (NEDZ - NBGZ) % NSKPZ;
    // number of recording points in total
    rec_NX = (NEDX - NBGX) / NSKPX + 1;
    rec_NY = (NEDY - NBGY) / NSKPY + 1;
    rec_NZ = (NEDZ - NBGZ) / NSKPZ + 1;

    // specific to each processor:
    calcRecordingPoints(&rec_nbgx, &rec_nedx, &rec_nbgy, &rec_nedy, &rec_nbgz, &rec_nedz,
            &rec_nxt, &rec_nyt, &rec_nzt, &displacement,
            (long int)ni, (long int)nj, (long int)nk, 
            rec_NX, rec_NY, rec_NZ, NBGX, NEDX, NSKPX, NBGY,
            NEDY, NSKPY, NBGZ, NEDZ, NSKPZ, thisid);

    int *disp = Alloc1P(size);
    int *rcount = Alloc1P(size);

    /// displacement, rec_nxt, rec_nyt, rec_nzt, WRITE_STEP
    unsigned int disp_xyzw_size = size * 5;
    unsigned int *disp_xyzw = Alloc1PU(disp_xyzw_size);

    unsigned int dxyzw[5] = {displacement, rec_nxt, rec_nyt, rec_nzt, WRITE_STEP};
    MPI_Allgather(dxyzw, 5, MPI_UNSIGNED, disp_xyzw, 5, MPI_UNSIGNED, MPI_COMM_WORLD);
    unsigned int max_send_size = 0;
    for (i = 0; i < size; i++) {
        unsigned int send_size = disp_xyzw[i*5+1] * disp_xyzw[i*5+2] * disp_xyzw[i*5+3] * disp_xyzw[i*5+4];
        max_send_size = max_send_size < send_size ? send_size : max_send_size;
    }
    float *snapshots = Alloc1D(max_send_size * size);

    if (this_rank == 0) LOG_INFO("Allocate buffers of #elements: %d\n", max_send_size);
    float *Bufx = Alloc1D(max_send_size);
    float *Bufy = Alloc1D(max_send_size);
    float *Bufz = Alloc1D(max_send_size);

    int RNX = ((NEDX - NBGX + NSKPX) / NSKPX);
    int RNY = ((NEDY - NBGY + NSKPY) / NSKPY);
    int RNZ = ((NEDZ - NBGZ + NSKPZ) / NSKPZ);
    if (this_rank == 0) write_metadata(OUT, "metadata", disp_xyzw, max_send_size,
            RNX, RNY, RNZ, WRITE_STEP, size * 5);

    //load checkpoint
     if(LOAD_CKPT) {
         load_ckpt(&it, PX, PY, PZ, this_rank, nx, ny, lnz, WSIZE, MSIZE, CUR_CKPT_IDX, 
                 NUM_CKPT_DIR, W_8, W_1, MPI_COMM_WORLD);
     }

    /* allocate buffer for buffering messages */
    void* buff_addr;
    int buffersize;
    buffersize = LenFDL * 8 * WSIZE * (_max(nx * lnz, _max(ny * lnz, nx * ny))) * sizeof(MPI_FLOAT);
    //MPI_Barrier(SWMPI_COMM);
    buff_addr = malloc(buffersize);
    if (!buff_addr) printf("allocation failure for buffer for MPI_Bsend !\n");
    MPI_Buffer_attach(buff_addr, buffersize);

    int load_ckpt_first = 1;
	float t, t0, t_weight;
	int srcIt = 0;
    while(1) {
        if(it >= NT) break;
        t_start = MPI_Wtime();
        // 8 steps as an one big loop
        int istep = 1;
        if(load_ckpt_first == 1 && LOAD_CKPT) {
            istep = it % 8 + 1;
            load_ckpt_first = 0;
        }
        for (; istep <=8; istep++){
            //Add_Source(W,M,it);
            ret = Gstart ("interpMomentRate");              // Start a manual timer

		    t = it * DT;
		    srcIt = int( t / src_in.dt );
	        if ( 0 != src_in.npts && srcIt < src_in.nt - 2 && params.useMultiSource )
	        {
		    	t0 = float( srcIt ) * src_in.dt;
		    	t_weight = ( t - t0 ) / src_in.dt;
                if (masternode)
                    printf( "srcIt = %d\n", srcIt );
		    	interpMomentRate( src_in, srcSlicePos, momentRate, momentRateSlice, t_weight, srcIt );
		    	// addMomentRate( W, grid, momentRateSlice, srcIndex, srcSlicePos, src_in.npts, DH, DT, M, it );
                addMomentRate(W_8, W_1, grid, momentRateSlice, srcIndex, srcSlicePos, src_in.npts, DH, DT, M, it );
		    }
            MPI_Barrier(MPI_COMM_WORLD);
            ret = Gstop ("interpMomentRate");              // Start a manual timer
            
            // mergeW(W, W_8, W_1);
            if (!params.useMultiSource) 
                AddSourceRicker(W_8,W_1,M,it);

            ret = Gstart ("RK_Syn");              // Start a manual timer
#ifdef ARCH_SW
            RK_Syn(W_8,W_1,M,hW_8,hW_1,mW_8,mW_1,tW_8,tW_1,matVx2Vz,matVy2Vz,istep,pml,isx1,isx2,isy1,isy2,isz1,isz2, &Aux);
#elif ARCH_x86
            RK_Syn(W,M,matVx2Vz,matVy2Vz,istep, pml, isx1, isx2, isy1, isy2, isz1, isz2, &Aux, DH, freenode);
#elif ARCH_MPE
            RK_Syn(W,M,matVx2Vz,matVy2Vz,istep, pml, isx1, isx2, isy1, isy2, isz1, isz2, &Aux, DH, freenode);
#endif
            ret = Gstop ("RK_Syn");              // Start a manual timer

            it += 1;
            /// save the sampled snapshot/*{{{*/
            int idtmp, tmpInd;
            if (it % TSKP == 0) {
                // added for plasticity
                idtmp = ((it / TSKP + WRITE_STEP - 1) % WRITE_STEP);
                idtmp = idtmp * rec_nxt * rec_nyt * rec_nzt;
                tmpInd = idtmp;
                // surface: k=nzt+align-1;
                int i, j, k;
                for (k = nz - ni1 - 1 - rec_nbgz; k >= nz - ni1 - 1 - rec_nedz; k = k - NSKPZ)
                    for (j = nj1 + rec_nbgy; j <= nj1 + rec_nedy; j = j + NSKPY)
                        for (i = ni1 + rec_nbgx; i <= ni1 + rec_nedx; i = i + NSKPX) {
                            // int idx = (i * ny * lnz + j * lnz + k) * WSIZE;
                            // Bufx[tmpInd] = W[idx + 0];
                            // Bufy[tmpInd] = W[idx + 1];
                            // Bufz[tmpInd] = W[idx + 2];

                            int idx = (i * ny * lnz + j * lnz + k) * WSIZE_V;
                            Bufx[tmpInd] = W_8[idx + 0];
                            Bufy[tmpInd] = W_8[idx + 1];
                            Bufz[tmpInd] = W_8[idx + 2];

                            tmpInd++;
                        }
               // save station info 
                if ( stationNum > 0 ) storageStation(grid, NT, stationNum, station, W_8, W_1, it);
                if (freenode == 1) comparePGV(grid, thisMPICoord, W_8, pgv);

                /// save snapshots
                if ((it / TSKP) % WRITE_STEP == 0) {
                    if (this_rank == 0) LOG_INFO("saving sampled snapshots @ timestep %ld", it);
                        // data2D_XYZ_out(thisMPICoord, params, grid, W_8, W_1, slice, sliceData, 'V', it, SAMP, SAMPLE_SIZE,sub_rank_x,sub_comm_x, sub_rank_y, sub_comm_y, Gather);
                    // write_model(this_rank, it, filenamebasex, Bufx, max_send_size, snapshots, max_send_size * size, MPI_COMM_WORLD);
                    // write_model(this_rank, it, filenamebasey, Bufy, max_send_size, snapshots, max_send_size * size, MPI_COMM_WORLD);
                    // write_model(this_rank, it, filenamebasez, Bufz, max_send_size, snapshots, max_send_size * size, MPI_COMM_WORLD);
                }
            } /// end of if (it % TSKP == 0)/*}}}*/

            //writeDat(W, it, this_rank);
            if(SAVE_FULL_IMG && it >= SAVE_FULL_IMG_START_STEP && (it - SAVE_FULL_IMG_START_STEP)  % SAVE_FULL_IMG_PER_STEP == 0) {
                int TARGET_OUTPUT_AT_DEPTH = 0; /// counting from 0
                write_full_snapshot(it, NX, NY, PX, PY, OUT, this_rank, MPI_COMM_WORLD, thisid, disp_xyzw, size, TARGET_OUTPUT_AT_DEPTH, W_8);
            }

            //save checkpoint
             if(SAVE_CKPT && it >= SAVE_CKPT_START_STEP && (it - SAVE_CKPT_START_STEP) % SAVE_CKPT_PER_STEP == 0) {
                 //save_ckpt(it, PX, PY, PZ, this_rank, nx, ny, lnz, WSIZE, MSIZE, CUR_CKPT_IDX, NUM_CKPT_DIR, W,hW, mW, tW, M, MPI_COMM_WORLD);
                 save_ckpt(it, PX, PY, PZ, this_rank, nx, ny, lnz, WSIZE, MSIZE, CUR_CKPT_IDX, NUM_CKPT_DIR, W_8, W_1, MPI_COMM_WORLD);
                 CUR_CKPT_IDX ++;
             }

            // if( it!=0 && it % 5000 == 0) export_seismo(R, rec_index, it);

        }// end of time 

        t_end = MPI_Wtime();
        if (masternode)
            printf("%6d of Total Timesteps = %6d;  cost %f sec  per step\n", it, NT, (t_end - t_start)/8.0f);

    }

	if ( stationNum > 0 ) {
		stationWrite(params, grid, thisMPICoord, station, NT, stationNum);
		freeStation(station);
    }

    if (freenode == 1) {
        outputPGV(params, grid, thisMPICoord, pgv);
		freePGV(pgv);
    }

    // export_seismo(R, rec_index, RNT);
    t_end = MPI_Wtime();
    if (masternode){
        int hr, mi;
        double totalT, sec;
        totalT = t_end - t_start0;
        hr = (int)(totalT / 3600.0);
        mi = (int)((totalT - hr*3600.0) / 60.0);
        sec = totalT - hr*3600.0 - mi*60.0;
        printf("Total time cost = %4dhour %2dmin %lfsec\n",
                hr, mi, sec);
    }
    if ( 0 != src_in.npts ) {
    	freeMomentRateSlice(momentRateSlice);
        free( srcIndex );
        free( srcSlicePos );
        freeReadMomentRate( momentRate );
    }

    /* de-allocate buffer for messages */
    MPI_Buffer_detach(buff_addr, &buffersize);
    free(wp_yzs0);
    free(wp_yzr0);
    free(wp_yzs1);
    free(wp_yzr1);
    free(wp_xzs0);
    free(wp_xzr0);
    free(wp_xzs1);
    free(wp_xzr1);
    free(wp_xys0);
    free(wp_xyr0);
    free(wp_xys1);
    free(wp_xyr1);

    free(buff_addr);
    free(C);
    //if(isource) { if(masternode) free(S); } 
    free(M);
    // free(W);
    // free(hW);
    // free(mW);
    // free(tW);
    free(W_8);
    free(hW_8);
    free(mW_8);
    free(tW_8);
    free(W_1);
    free(hW_1);
    free(mW_1);
    free(tW_1);
    free(rec_index);
    free(R);
    Delloc3D(damp);   

    ret = Gstop("main");
    // ret = Gstop ("main");
    // char path[256];
    // sprintf(path, "./timing/%s", this_rank);
    // ret = GPTLpr_file(path);
    GPTLpr(this_rank);
    if (GPTLfinalize () != 0)
       return 1;

    // MPI_Finalize();
    return 0;
}
