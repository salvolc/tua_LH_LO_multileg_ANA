/*
VectorXd vertex_match(TClonesArray* TCP){
//	int npart_sel = 50;
	int nPart = TCP->GetEntries();
	MatrixXd particle_infos(nPart,7);

	int nt = 0;int i_nt = 0;
	int nu = 0;int i_nu = 0;
	int ng = 0;int i_ng = 0;
	int ngam = 0;int i_ngam = 0;

	for (int iPart = 0; iPart < nPart; ++iPart)
	{

		GenParticle* P_Particle = (GenParticle*)TCP->At(iPart);
		GenParticle* P_ParticleHelp = (GenParticle*)TCP->At(iPart);

		particle_infos(iPart,0) = P_Particle->PID;
		if(abs(P_Particle->M1) == 1){particle_infos(iPart,1) = 0;}
		else{
			P_ParticleHelp = (GenParticle*)TCP->At(abs(P_Particle->M1));
			particle_infos(iPart,1) = P_ParticleHelp->PID;
		}
		if(abs(P_Particle->M2) == 1){particle_infos(iPart,2) = 0;}
		else{
			P_ParticleHelp = (GenParticle*)TCP->At(abs(P_Particle->M2));
			particle_infos(iPart,2) = P_ParticleHelp->PID;
		}
		if(abs(P_Particle->D1) == 1){particle_infos(iPart,3) = 0;}
		else{
			P_ParticleHelp = (GenParticle*)TCP->At(abs(P_Particle->D1));
			particle_infos(iPart,3) = P_ParticleHelp->PID;
		}
		if(abs(P_Particle->D2) == 1){particle_infos(iPart,4) = 0;}
		else{
			P_ParticleHelp = (GenParticle*)TCP->At(abs(P_Particle->D2));
			particle_infos(iPart,4) = P_ParticleHelp->PID;
		}
		particle_infos(iPart,5) = iPart;
		particle_infos(iPart,6) = abs(P_Particle->PT);
		if(abs(particle_infos(iPart,0)) == 6)nt++;
		if(abs(particle_infos(iPart,0)) == 2 || abs(particle_infos(iPart,0)) == 4)nu++;
		if(abs(particle_infos(iPart,0)) == 21 || abs(particle_infos(iPart,0)) == 9)ng++;
		if(abs(particle_infos(iPart,0)) == 22)ngam++;
	}


	MatrixXd t_infos = MatrixXd::Zero(nt,7);
	MatrixXd u_infos = MatrixXd::Zero(nu,7);
	MatrixXd g_infos = MatrixXd::Zero(ng,7);
	MatrixXd gam_infos = MatrixXd::Zero(ngam,7);

	for (int iPart = 0; iPart < nPart; ++iPart)
	{
		if(abs(particle_infos(iPart,0)) == 6){t_infos.row(i_nt) = particle_infos.row(iPart);i_nt++;}
		if(abs(particle_infos(iPart,0)) == 22){gam_infos.row(i_ngam) = particle_infos.row(iPart);i_ngam++;}
		if(abs(particle_infos(iPart,0)) == 2 || abs(particle_infos(iPart,0)) == 4){u_infos.row(i_nu) = particle_infos.row(iPart);i_nu++;}
		if(abs(particle_infos(iPart,0)) == 21 || abs(particle_infos(iPart,0)) == 9){g_infos.row(i_ng) = particle_infos.row(iPart);i_ng++;}
	}

	if(debug)cout << particle_infos.block(0,0,100,7) << endl;

	VectorXd i_vertex_matched = VectorXd::Zero(8);
	bool matched_t = false; 
	bool matched_u = false; i_nu=i_nu/2;
	bool matched_b = false; i_ngam=i_ngam/2;
	bool matched_g = false; i_ng=i_ng/2; 
	GenParticle* P_ParticleHelp = (GenParticle*)TCP->At(0);
	for (int it = 0; it < i_nt && !matched_t; ++it)
	{
		if ( (t_infos(it,3)==22) && (abs(t_infos(it,4))==4  ||  abs(t_infos(it,4))==2 ) )   
		{
			i_vertex_matched(1) = t_infos(it,5);
			P_ParticleHelp = (GenParticle*)TCP->At(i_vertex_matched(1));
			i_vertex_matched(0) = P_ParticleHelp->D1;
			i_vertex_matched(2) = P_ParticleHelp->D2;
			matched_t = true;
		}
		if ( (t_infos(it,4)==22) && (abs(t_infos(it,3))==4  ||  abs(t_infos(it,3))==2 ) )
		{
			i_vertex_matched(1) = t_infos(it,5);
			P_ParticleHelp = (GenParticle*)TCP->At(i_vertex_matched(1));
			i_vertex_matched(0) = P_ParticleHelp->D2;
			i_vertex_matched(2) = P_ParticleHelp->D1;
			matched_t = true;
		}
		if ( (abs(t_infos(it,4))==24) && (abs(t_infos(it,3))==5))
		{
			i_vertex_matched(3) = t_infos(it,5);
			P_ParticleHelp = (GenParticle*)TCP->At(i_vertex_matched(3));
			i_vertex_matched(4) = P_ParticleHelp->D1;
			i_vertex_matched(5) = P_ParticleHelp->D2;
			matched_b = true;
		}
		if ( (abs(t_infos(it,3))==24) && (abs(t_infos(it,4))==5))
		{
			i_vertex_matched(3) = t_infos(it,5);
			P_ParticleHelp = (GenParticle*)TCP->At(i_vertex_matched(3));
			i_vertex_matched(4) = P_ParticleHelp->D2;
			i_vertex_matched(5) = P_ParticleHelp->D1;
			matched_b = true;
		}
	}

	if(!matched_t){//Falls t nicht so gefunden wird haben wir wohl ein up was nach t gamma geht
	for (int iu = 0; iu < i_nu && !matched_u; ++iu)
	{
		if ( (u_infos(iu,3)==22) && (abs(u_infos(iu,4))==6) )   
		{
			i_vertex_matched(2) = u_infos(iu,5);
			P_ParticleHelp = (GenParticle*)TCP->At(i_vertex_matched(2));
			i_vertex_matched(0) = P_ParticleHelp->D1;
			i_vertex_matched(1) = P_ParticleHelp->D2;
			matched_u = true;
		}
		if ( (u_infos(iu,4)==22) && (abs(u_infos(iu,3))==6) )
		{
			i_vertex_matched(2) = u_infos(iu,5);
			P_ParticleHelp = (GenParticle*)TCP->At(i_vertex_matched(2));
			i_vertex_matched(0) = P_ParticleHelp->D2;
			i_vertex_matched(1) = P_ParticleHelp->D1;
			matched_u = true;
		}
	}
	}
	if(!matched_t && !matched_u){
	for (int ig = 0; ig < i_ng && !matched_g; ++ig)
	{
		if ( (abs(g_infos(ig,3))==6) && ( (abs(g_infos(ig,4))==4) || (abs(g_infos(ig,4))==2) ) )   
		{
			i_vertex_matched(7) = g_infos(ig,5);
			P_ParticleHelp = (GenParticle*)TCP->At(i_vertex_matched(7));
			i_vertex_matched(1) = P_ParticleHelp->D1;
			i_vertex_matched(2) = P_ParticleHelp->D2;
			matched_g = true;
		}
		if ( abs((g_infos(ig,4))==6) && ( (abs(g_infos(ig,3))==4) || (abs(g_infos(ig,3))==2) ) )
		{
			i_vertex_matched(7) = g_infos(ig,5);
			P_ParticleHelp = (GenParticle*)TCP->At(i_vertex_matched(7));
			i_vertex_matched(1) = P_ParticleHelp->D2;
			i_vertex_matched(2) = P_ParticleHelp->D1;
			matched_g = true;
		}
	}
	for (int igam = 0; igam < i_ngam; ++igam)
	{
		if ( (abs(gam_infos(igam,1))==21  ||  abs(gam_infos(igam,1))==9 ) || (abs(gam_infos(igam,2))==21  ||  abs(gam_infos(igam,2))==9 ) )   
		{
			i_vertex_matched(0) = gam_infos(igam,5);
		}
	}
	}
	if(matched_u && debug)cout << "Up gemached" << endl;
	if(matched_t && debug)cout << "Top gemached" << endl;
	if(matched_g && debug)cout << "Gluon gemached" << endl;
	if(matched_t || matched_u || matched_g){i_vertex_matched(6)=1;}else{i_vertex_matched(6)=0;}
	if(matched_b){i_vertex_matched(7)=1;}else{i_vertex_matched(7)=0;}
	if(debug){cout << "Der Vertex Vektor:"<< endl;cout << i_vertex_matched << endl;}
	return i_vertex_matched;
}
*/