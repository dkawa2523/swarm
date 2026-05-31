# Fixed-particle MC warmup comparison

Warmup: 600 collisions. Production sampling: 1200 collisions. Particles: 768.

case_id,E_over_N_Td,mean_energy_eV,drift_velocity_m_s,net_ionization_frequency_s,mc_population_model,ionization_source_treatment,ionization_branching_model,secondary_electron_tracking,mc_population_control,mc_resampling_events,mc_effective_particle_count,inelastic_angular_model,mc_energy_balance_status,mc_tracked_energy_balance_residual_fraction,mc_weight_balance_residual_fraction,mc_physical_branching_gap_eV,mc_null_collision_acceptance_fraction,mc_max_collision_to_trial_ratio,mc_tail_uncertainty_status,mc_min_tail_bin_count,mc_tail_effective_sample_count_min,mc_tail_weak_probability_fraction,mc_max_resolved_energy_eV,mc_tail_comparison_status
fixed_warmup_hp_0000,100.0,6.503388155601828,226126.23732484065,542634.7867443362,fixed_particle_single_daughter,equal,single_daughter_sampling,False,none,0,768.0,isotropic_reset,ok,0.000194950324159,0.0,221.94081018067163,0.040884693287037,0.1411856837004381,ok,5.746145514685116,5.746145514685116,0.0116694497616712,23.25,ok
fixed_warmup_hp_0001,200.0,8.012756697265818,380105.72451935057,4161695.567664051,fixed_particle_single_daughter,equal,single_daughter_sampling,False,none,0,768.0,isotropic_reset,ok,8.139604074236229e-05,0.0,3117.603581049366,0.0492607060185185,0.187323828032005,ok,1.867200385737978,1.867200385737978,0.0041092917998614,33.25,ok
fixed_warmup_hp_0002,400.0,10.813234094267257,603910.0353605595,18188906.416874997,fixed_particle_single_daughter,equal,single_daughter_sampling,False,none,0,768.0,isotropic_reset,ok,4.007349381847617e-05,0.0,26415.832357159496,0.0621925636574074,0.2402791261655344,ok,1.0,1.0,0.0027178759474811,50.75,ok
fixed_warmup_hp_0003,600.0,13.742013988033426,773950.3553129677,37398634.73073806,fixed_particle_single_daughter,equal,single_daughter_sampling,False,none,0,768.0,isotropic_reset,ok,2.899308599320733e-05,0.0,78083.70537636586,0.0740357349537037,0.2815648705106494,ok,5.045651214639339,5.045651214639339,0.0022078609777176,105.0,ok

## Metrics

E_over_N_Td,reference,candidate,eedf_relative_l1,reference_mean_energy_eV,candidate_mean_energy_eV,mean_energy_relative_difference,mc_energy_balance_status,mc_tail_comparison_status,mc_tail_weak_probability_fraction,mc_max_resolved_energy_eV
100.0,two_term native SG hold,MCIG GUI actual,0.017114680313148476,6.496528329508164,6.427562525573979,0.010615793611017258,ok,ok,0.0116694497616712,23.25
100.0,two_term native SG hold,"internal MC fixed warmup HP (768, warmup 600 + production 1200)",0.023595097112309938,6.496528329508164,6.546843971347718,0.007745004606693166,ok,ok,0.0116694497616712,23.25
100.0,MCIG GUI actual,"internal MC fixed warmup HP (768, warmup 600 + production 1200)",0.03142646655097283,6.427562525573979,6.546843971347718,0.018557804035222045,ok,ok,0.0116694497616712,23.25
200.0,two_term native SG hold,MCIG GUI actual,0.04205761611603273,7.953861753381179,7.757728733922895,0.02465884189839075,ok,ok,0.0041092917998614,33.25
200.0,two_term native SG hold,"internal MC fixed warmup HP (768, warmup 600 + production 1200)",0.014072044368871174,7.953861753381179,8.062217393895653,0.013623022862876892,ok,ok,0.0041092917998614,33.25
200.0,MCIG GUI actual,"internal MC fixed warmup HP (768, warmup 600 + production 1200)",0.04526834256512883,7.757728733922895,8.062217393895653,0.0392497173356028,ok,ok,0.0041092917998614,33.25
400.0,two_term native SG hold,MCIG GUI actual,0.082517932641006,10.56937493055106,10.10913280279671,0.0435448766628581,ok,ok,0.0027178759474811,50.75
400.0,two_term native SG hold,"internal MC fixed warmup HP (768, warmup 600 + production 1200)",0.026179314492431317,10.56937493055106,10.85694522113558,0.02720788054866811,ok,ok,0.0027178759474811,50.75
400.0,MCIG GUI actual,"internal MC fixed warmup HP (768, warmup 600 + production 1200)",0.08540709908910667,10.10913280279671,10.85694522113558,0.07397394345556398,ok,ok,0.0027178759474811,50.75
600.0,two_term native SG hold,MCIG GUI actual,0.1103998407733816,13.163543564923732,12.461445029668953,0.05333659069778347,ok,ok,0.0022078609777176,105.0
600.0,two_term native SG hold,"internal MC fixed warmup HP (768, warmup 600 + production 1200)",0.05056944532456399,13.163543564923732,13.782177622115976,0.04699601244460411,ok,ok,0.0022078609777176,105.0
600.0,MCIG GUI actual,"internal MC fixed warmup HP (768, warmup 600 + production 1200)",0.10799700297074484,12.461445029668953,13.782177622115976,0.10598550884769331,ok,ok,0.0022078609777176,105.0
