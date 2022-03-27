#!/usr/bin/env python

"""External libaries"""
import pickle
import analysis

# Results_0 = pickle.load(open( "Results_0.p", "rb" ) )
# Results_1 = pickle.load(open( "Results_1.p", "rb" ) )
# Results_2 = pickle.load(open( "Results_2.p", "rb" ) )
# Results = analysis.concatenate_results(Results_0, Results_1, Results_2)

Results = pickle.load(open( "Results_0.p", "rb" ) )

print(analysis.calculate_mean_drift_velocity(Results))

analysis.plot_results(Results)
