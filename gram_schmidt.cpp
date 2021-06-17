void GramSchmidt(){
	vector<vector<double>> v;
	v.push_back({5, 6, 3});
	v.push_back({3, 2, 6});
	v.push_back({4, 5, 8});
	double norm_v1 = norm(v[0]);
	for(int i = 0; i < v[0].size(); i++){
		v[0][i] = v[0][i]/norm_v1;
	}
	for(int i = 1; i < v.size(); i++){
		vector<double> w = v[i];
		for(int j = 0; j < i; j++){
			double proj = inner_product(w.begin(), w.end(), v[j].begin(), 0.0)/1;
			for(int k = 0; k < v[i].size(); k++){
				w[k] = w[k] - proj*v[j][k];
			}
		}
		double norm_w = norm(w);
		for(int l = 0; l < w.size(); l++){
			v[i][l] = w[l]/norm_w;	
		}
	}

	for(int i = 0; i < v.size(); i++){
		printVector("v[i]", v[i]);
	}
	
}