%% Producing the main figures in Sheintuch et al., 2023:
% Illustrations are not shown here.
clc
clear all
close all

% Defining the parameters:
maximal_number_of_sessions=16;
num_trials=2;
days_vec=[1:2:15 19:2:33];
elapsed_days_vec=0:2:30;
novel_sessions=[1:3,9:11];
bin_size=4;
num_bins=20;

data_pathway='E:\downloads\cell_assemblies-main\Data\';
CA1_group={'C5M1','C6M3','C6M4','C8M2'};
CA3_group={'C14M4','C15M2','C23M4','C24M3','C24M4'};

% Loading the data:
list=dir(data_pathway);
a=1:numel(list);
names_of_files = arrayfun(@(x) list(x).name,a,'UniformOutput',false);
miceFilesIdx=find(cellfun (@(x) ~isempty(strfind(x,'C')) ,names_of_files));
number_of_mice=numel(miceFilesIdx);
across_mice_data=cell(1,number_of_mice);
mouse_group=zeros(1,number_of_mice);
CA1_ordered_names=cell(size(CA1_group));
CA3_ordered_names=cell(size(CA3_group));
CA1_count=0;
CA3_count=0;
for n=1:number_of_mice
    temp_data=load([data_pathway names_of_files{miceFilesIdx(n)}]);
    across_mice_data{n}=temp_data.this_mouse_data;
    
%     this_mouse_data=struct;
%     this_mouse_data.SI_SSR_bit_spike_sessions=across_mice_data{n}.SI_SSR_bit_spike_sessions;
%     this_mouse_data.average_SI_SSR_bit_spike_sessions=across_mice_data{n}.average_SI_SSR_bit_spike_sessions;
%     this_mouse_data.PV_correlations_bin_resolution=across_mice_data{n}.PV_correlations_bin_resolution;
%     this_mouse_data.within_day_PV_correlations=across_mice_data{n}.within_day_PV_correlations;
%     this_mouse_data.mean_PV_dynamics=across_mice_data{n}.mean_PV_dynamics;
%     this_mouse_data.mean_rate_dynamics=across_mice_data{n}.mean_rate_dynamics;
%     this_mouse_data.rate_correlations=across_mice_data{n}.rate_correlations;
%     this_mouse_data.PV_correlations=across_mice_data{n}.PV_correlations;
%     this_mouse_data.maximal_PV_correlations=across_mice_data{n}.maximal_PV_correlations;
%     this_mouse_data.number_of_neurons_sessions=across_mice_data{n}.number_of_neurons_sessions;
%     this_mouse_data.average_firing_rates_sessions=across_mice_data{n}.average_firing_rates_sessions;
%     this_mouse_data.number_of_track_traversals_sessions=across_mice_data{n}.number_of_track_traversals_sessions;
%     this_mouse_data.all_within_day_pairwise_correlations_per_mouse_novel=across_mice_data{n}.all_within_day_pairwise_correlations_per_mouse_novel;
%     this_mouse_data.same_bin_within_day_pairwise_correlations_per_mouse_novel=across_mice_data{n}.same_bin_within_day_pairwise_correlations_per_mouse_novel;
%     this_mouse_data.same_bin_across_sides_pairwise_correlations_per_mouse_novel=across_mice_data{n}.same_bin_across_sides_pairwise_correlations_per_mouse_novel;
%     this_mouse_data.diff_bins_within_day_pairwise_correlations_per_mouse_novel=across_mice_data{n}.diff_bins_within_day_pairwise_correlations_per_mouse_novel;
%     this_mouse_data.diff_bins_across_sides_pairwise_correlations_per_mouse_novel=across_mice_data{n}.diff_bins_across_sides_pairwise_correlations_per_mouse_novel;
%     
%     this_mouse_data.average_SI_SSR_bit_spike_dynamics_correlated_cells=across_mice_data{n}.average_SI_SSR_bit_spike_dynamics_correlated_cells;
%     this_mouse_data.average_SI_SSR_bit_spike_dynamics_uncorrelated_cells=across_mice_data{n}.average_SI_SSR_bit_spike_dynamics_uncorrelated_cells;
%     this_mouse_data.PV_correlations_correlated_cells=across_mice_data{n}.PV_correlations_correlated_cells;
%     this_mouse_data.PV_correlations_uncorrelated_cells_control=across_mice_data{n}.PV_correlations_uncorrelated_cells_control;
%     this_mouse_data.PV_dynamics_correlated_cells=across_mice_data{n}.PV_dynamics_correlated_cells;
%     this_mouse_data.PV_dynamics_uncorrelated_cells_control=across_mice_data{n}.PV_dynamics_uncorrelated_cells_control;
%     
%     this_mouse_data.within_day_tuning_corr_single_cells=across_mice_data{n}.within_day_tuning_corr_single_cells;
%     this_mouse_data.SI_SSR_single_cells=across_mice_data{n}.SI_SSR_single_cells;
%     this_mouse_data.correlation_percentile_single_cells=across_mice_data{n}.correlation_percentile_single_cells;
%     this_mouse_data.across_days_tuning_corr_single_cells=across_mice_data{n}.across_days_tuning_corr_single_cells;
%     this_mouse_data.across_days_correlation_percentile_single_cells=across_mice_data{n}.across_days_correlation_percentile_single_cells;
%     this_mouse_data.across_days_SI_SSR_single_cells=across_mice_data{n}.across_days_SI_SSR_single_cells;
%     
%     this_mouse_data.significant_right_trials=across_mice_data{n}.significant_right_trials;
%     this_mouse_data.significant_left_trials=across_mice_data{n}.significant_left_trials;
%     this_mouse_data.rate_maps_trials=across_mice_data{n}.rate_maps_trials;
%     
%     this_mouse_data.SI_SSR_bit_spike_sessions_subsampled=across_mice_data{n}.SI_SSR_bit_spike_sessions_subsampled;
%     this_mouse_data.average_SI_SSR_bit_spike_sessions_subsampled=across_mice_data{n}.average_SI_SSR_bit_spike_sessions_subsampled;
%     this_mouse_data.PV_correlations_bin_resolution_subsampled=across_mice_data{n}.PV_correlations_bin_resolution_subsampled;
%     this_mouse_data.within_day_PV_correlations_subsampled=across_mice_data{n}.within_day_PV_correlations_subsampled;
%     if n==1 || n==8
%         this_mouse_data.example_day_mouse_position=across_mice_data{n}.example_day_mouse_position;
%         this_mouse_data.example_day_prior_distribution=across_mice_data{n}.example_day_prior_distribution;
%         this_mouse_data.example_day_mouse_velocities=across_mice_data{n}.example_day_mouse_velocities;
%         this_mouse_data.example_cells_spiking_activity=across_mice_data{n}.example_cells_spiking_activity;
%         this_mouse_data.example_day_tuning_curves=across_mice_data{n}.example_day_tuning_curves;
%         this_mouse_data.example_day_concatanation_index=across_mice_data{n}.example_day_concatanation_index;
%     end
%     this_mouse_data.number_of_neurons_dynamics=across_mice_data{n}.number_of_neurons_dynamics;
%     this_mouse_data.cell_to_index_map=across_mice_data{n}.cell_to_index_map;
%     this_mouse_data.overall_registration_errors=across_mice_data{n}.overall_registration_errors;
%     this_mouse_data.PV_dynamics_subsampled=across_mice_data{n}.PV_dynamics_subsampled;
%     this_mouse_data.rate_dynamics_subsampled=across_mice_data{n}.rate_dynamics_subsampled;
%     
%     this_mouse_data.average_rates_correlated_sessions=across_mice_data{n}.average_rates_correlated_sessions;
%     this_mouse_data.average_rates_uncorrelated_sessions=across_mice_data{n}.average_rates_uncorrelated_sessions;
%     this_mouse_data.rate_correlations_correlated_cells=across_mice_data{n}.rate_correlations_correlated_cells;
%     this_mouse_data.rate_correlations_uncorrelated_cells=across_mice_data{n}.rate_correlations_uncorrelated_cells;
%     this_mouse_data.noise_correlations_same_pos_correlated_single_cells=across_mice_data{n}.noise_correlations_same_pos_correlated_single_cells;
%     this_mouse_data.noise_correlations_same_pos_uncorrelated_single_cells=across_mice_data{n}.noise_correlations_same_pos_uncorrelated_single_cells;
%     this_mouse_data.all_noise_correlations_same_pos_corr_cells_per_mouse=across_mice_data{n}.noise_correlations_same_pos_correlated_single_cells;
%     this_mouse_data.all_noise_correlations_same_pos_uncorr_cells_per_mouse=across_mice_data{n}.noise_correlations_same_pos_uncorrelated_single_cells;
%     
%     this_mouse_data.mean_significant_cells_is_sig_dynamics=across_mice_data{n}.mean_significant_cells_is_sig_dynamics;
%     this_mouse_data.mean_not_significant_cells_is_sig_dynamics=across_mice_data{n}.mean_not_significant_cells_is_sig_dynamics;
%     this_mouse_data.mean_significant_cells_is_sig_across_env=across_mice_data{n}.mean_significant_cells_is_sig_across_env;
%     this_mouse_data.mean_not_significant_cells_is_sig_across_env=across_mice_data{n}.mean_not_significant_cells_is_sig_across_env;
%     this_mouse_data.within_day_pairwise_tuning_corr_single_cells=across_mice_data{n}.within_day_pairwise_tuning_corr_single_cells;
%     this_mouse_data.within_day_SI_SSR_single_cells=across_mice_data{n}.within_day_SI_SSR_single_cells;
%     this_mouse_data.within_day_within_day_tuning_corr_single_cells=across_mice_data{n}.within_day_within_day_tuning_corr_single_cells;
%     this_mouse_data.acorss_days_across_days_tuning_corr_single_cells=across_mice_data{n}.acorss_days_across_days_tuning_corr_single_cells;
%     this_mouse_data.across_days_average_rates_single_cells=across_mice_data{n}.across_days_average_rates_single_cells;
%     this_mouse_data.within_day_pairwise_tuning_corr_single_cells_separated=across_mice_data{n}.within_day_pairwise_tuning_corr_single_cells_separated;
%     this_mouse_data.SI_SSR_single_cells_separated=across_mice_data{n}.SI_SSR_single_cells_separated;
%     this_mouse_data.within_day_average_rates_single_cells=across_mice_data{n}.within_day_average_rates_single_cells;
%     
%     this_mouse_data.all_pairwise_correlations_session_one_right_across=across_mice_data{n}.all_pairwise_correlations_session_one_right_across;
%     this_mouse_data.all_pairwise_correlations_session_one_left_across=across_mice_data{n}.all_pairwise_correlations_session_one_left_across;
%     this_mouse_data.all_pairwise_correlations_session_one_right_within=across_mice_data{n}.all_pairwise_correlations_session_one_right_within;
%     this_mouse_data.all_pairwise_correlations_session_one_left_within=across_mice_data{n}.all_pairwise_correlations_session_one_left_within;
%     this_mouse_data.all_pairwise_correlations_session_two_right_across=across_mice_data{n}.all_pairwise_correlations_session_two_right_across;
%     this_mouse_data.all_pairwise_correlations_session_two_left_across=across_mice_data{n}.all_pairwise_correlations_session_two_left_across;
%     this_mouse_data.all_pairwise_correlations_session_two_right_within=across_mice_data{n}.all_pairwise_correlations_session_two_right_within;
%     this_mouse_data.all_pairwise_correlations_session_two_left_within=across_mice_data{n}.all_pairwise_correlations_session_two_left_within;
%     this_mouse_data.average_noise_corr_versus_shift_novel=across_mice_data{n}.average_noise_corr_versus_shift_novel;
%     this_mouse_data.shuffle_noise_correlation_versus_shift_novel=across_mice_data{n}.shuffle_noise_correlation_versus_shift_novel;
%     this_mouse_data.average_noise_correlation_versus_distance_novel=across_mice_data{n}.average_noise_correlation_versus_distance_novel;
%     this_mouse_data.average_shuffle_noise_correlation_versus_distance_novel=across_mice_data{n}.average_shuffle_noise_correlation_versus_distance_novel;
%     this_mouse_data.maximal_PV_correlations=across_mice_data{n}.maximal_PV_correlations;
%     
%     this_mouse_data.average_number_of_fields_dynamics_dependent=across_mice_data{n}.average_number_of_fields_dynamics_dependent;
%     this_mouse_data.average_max_field_size_dynamics_dependent=across_mice_data{n}.average_max_field_size_dynamics_dependent;
%     this_mouse_data.average_number_of_fields_dynamics_independent=across_mice_data{n}.average_number_of_fields_dynamics_independent;
%     this_mouse_data.average_max_field_size_dynamics_independent=across_mice_data{n}.average_max_field_size_dynamics_independent;
%     
%     this_mouse_data.same_bin_within_day_pairwise_correlations_per_mouse_novel_CTC=across_mice_data{n}.same_bin_within_day_pairwise_correlations_per_mouse_novel_CTC;
%     this_mouse_data.same_bin_across_mice_pairwise_correlations_per_mouse_novel_CTC=across_mice_data{n}.same_bin_across_mice_pairwise_correlations_per_mouse_novel_CTC;
%     this_mouse_data.same_bin_across_sides_pairwise_correlations_per_mouse_novel_CTC=across_mice_data{n}.same_bin_across_sides_pairwise_correlations_per_mouse_novel_CTC;
%     this_mouse_data.same_bin_within_day_pairwise_correlations_per_mouse_novel_AM=across_mice_data{n}.same_bin_within_day_pairwise_correlations_per_mouse_novel_AM;
%     this_mouse_data.same_bin_across_mice_pairwise_correlations_per_mouse_novel_AM=across_mice_data{n}.same_bin_across_mice_pairwise_correlations_per_mouse_novel_AM;
%     
%     this_mouse_data.all_correlated_preferred_positions_right=all_correlated_preferred_positions_right{n};
%     this_mouse_data.all_correlated_preferred_positions_left=all_correlated_preferred_positions_left{n};
%     this_mouse_data.all_uncorrelated_preferred_positions_right=all_uncorrelated_preferred_positions_right{n};
%     this_mouse_data.all_uncorrelated_preferred_positions_left=all_uncorrelated_preferred_positions_left{n};
%     
    this_group=[];
    this_mouse=names_of_files{miceFilesIdx(n)};
    this_mouse=erase(this_mouse,'.mat');
    across_mice_data{n}.mouse_name=this_mouse;
    if find(cellfun (@(x) ~isempty(strfind(x,this_mouse)) ,CA1_group))
        this_group=1;
        CA1_count=CA1_count+1;
        CA1_ordered_names{CA1_count}=this_mouse;
    elseif find(cellfun (@(x) ~isempty(strfind(x,this_mouse)) ,CA3_group))
        this_group=2;
        CA3_count=CA3_count+1;
        CA3_ordered_names{CA3_count}=this_mouse;
    end
    across_mice_data{n}.group=this_group;
    mouse_group(n)=this_group;
    
%     save(this_mouse,'this_mouse_data');
end
CA1_indexes=find(mouse_group==1);
CA3_indexes=find(mouse_group==2);

%% Figure 1 - Example place cells imaged in hippocampal CA1 and CA3 during linear track familiarization:

% Figure 1C-D - Example place cells:
exmple_mouse_CA1='C6M4';
exmple_mouse_CA3='C14M4';
fs=20;
total_num_bins=24;
for n=1:number_of_mice
    if strcmp(across_mice_data{n}.mouse_name,exmple_mouse_CA1)
        exmple_mouse_index_CA1=n;
    elseif strcmp(across_mice_data{n}.mouse_name,exmple_mouse_CA3)
        exmple_mouse_index_CA3=n;
    end
end

% Loading the necessary data:
example_day_mouse_position_CA1=across_mice_data{exmple_mouse_index_CA1}.example_day_mouse_position;
example_day_prior_distribution_CA1=across_mice_data{exmple_mouse_index_CA1}.example_day_prior_distribution;
example_day_mouse_velocities_CA1=across_mice_data{exmple_mouse_index_CA1}.example_day_mouse_velocities;
example_cells_spiking_activity_CA1=across_mice_data{exmple_mouse_index_CA1}.example_cells_spiking_activity;
example_day_tuning_curves_CA1=across_mice_data{exmple_mouse_index_CA1}.example_day_tuning_curves;
example_day_concatanation_index_CA1=across_mice_data{exmple_mouse_index_CA1}.example_day_concatanation_index;
t_vec_CA1=0:1/fs:(length(example_day_mouse_position_CA1)-1)/fs;

example_day_mouse_position_CA3=across_mice_data{exmple_mouse_index_CA3}.example_day_mouse_position;
example_day_prior_distribution_CA3=across_mice_data{exmple_mouse_index_CA3}.example_day_prior_distribution;
example_day_mouse_velocities_CA3=across_mice_data{exmple_mouse_index_CA3}.example_day_mouse_velocities;
example_cells_spiking_activity_CA3=across_mice_data{exmple_mouse_index_CA3}.example_cells_spiking_activity;
example_day_tuning_curves_CA3=across_mice_data{exmple_mouse_index_CA3}.example_day_tuning_curves;
example_day_concatanation_index_CA3=across_mice_data{exmple_mouse_index_CA3}.example_day_concatanation_index;
t_vec_CA3=0:1/fs:(length(example_day_mouse_position_CA3)-1)/fs;

% Figure 1C - Example place cells in CA1:
figure('units','normalized','outerposition',[0.2 0.1 0.6 0.8])
for n=1:size(example_cells_spiking_activity_CA1,2)
    axes('position',[(n-1)*0.185+0.07 0.65 0.165 0.24])
    plot(bin_size*example_day_mouse_position_CA1,t_vec_CA1(end)-t_vec_CA1)
    chosen_cell_spikes=find(example_cells_spiking_activity_CA1(:,n)>0);
    chosen_cell_velocities=example_day_mouse_velocities_CA1(chosen_cell_spikes);
    
    num_spikes_chosen_cell=length(chosen_cell_spikes);
    for k=1:num_spikes_chosen_cell
        hold on
        if chosen_cell_velocities(k)>0
            plot(bin_size*example_day_mouse_position_CA1(chosen_cell_spikes(k)),t_vec_CA1(example_day_concatanation_index_CA1)-t_vec_CA1(chosen_cell_spikes(k)),'.','markersize',20,'color','r')
        else
            plot(bin_size*example_day_mouse_position_CA1(chosen_cell_spikes(k)),t_vec_CA1(example_day_concatanation_index_CA1)-t_vec_CA1(chosen_cell_spikes(k)),'.','markersize',20,'color',[0 0.8 0])
        end
        ylim([0 t_vec_CA1(example_day_concatanation_index_CA1)])
    end
    if n==1
        ylabel('Time (sec)')
    end
    set(gca,'ytick',[])
    set(gca,'xtick',[])
    hold on
    title(['Cell ' num2str(n)],'fontweight','normal')
    set(gca,'fontsize',16)
    if n==3
        text(50,800,'Figure 1C - CA1','fontsize',20,'fontweight','bold','HorizontalAlignment','Center')
    end
    
    axes('position',[(n-1)*0.185+0.07 0.39 0.165 0.24])
    plot(bin_size*example_day_mouse_position_CA1(example_day_concatanation_index_CA1+1:end),t_vec_CA1(example_day_concatanation_index_CA1)+t_vec_CA1(end)-t_vec_CA1(example_day_concatanation_index_CA1+1:end))
    chosen_cell_velocities=example_day_mouse_velocities_CA1(chosen_cell_spikes);
    
    num_spikes_chosen_cell=length(chosen_cell_spikes);
    for k=1:num_spikes_chosen_cell
        hold on
        if chosen_cell_velocities(k)>0
            plot(bin_size*example_day_mouse_position_CA1(chosen_cell_spikes(k)),t_vec_CA1(example_day_concatanation_index_CA1)+t_vec_CA1(end)-t_vec_CA1(chosen_cell_spikes(k)),'.','markersize',20,'color','r')
        else
            plot(bin_size*example_day_mouse_position_CA1(chosen_cell_spikes(k)),t_vec_CA1(example_day_concatanation_index_CA1)+t_vec_CA1(end)-t_vec_CA1(chosen_cell_spikes(k)),'.','markersize',20,'color',[0 0.8 0])
        end
        ylim([t_vec_CA1(example_day_concatanation_index_CA1+1) t_vec_CA1(end)])
    end
    if n==1
        ylabel('Time (sec)')
    end
    set(gca,'ytick',[])
    set(gca,'xtick',[])
    hold on
    set(gca,'fontsize',16)
    
    axes('position',[(n-1)*0.185+0.07 0.08 0.165 0.27])
    plot(-bin_size/2+bin_size*[1:total_num_bins],fs*(squeeze(example_day_tuning_curves_CA1(n,1,:)))/length(example_day_mouse_position_CA1)./squeeze(example_day_prior_distribution_CA1(:,1)),'color','r','linewidth',2)
    hold on
    plot(-bin_size/2+bin_size*[1:total_num_bins],fs*(squeeze(example_day_tuning_curves_CA1(n,2,:)))/length(example_day_mouse_position_CA1)./squeeze(example_day_prior_distribution_CA1(:,2)),'color','g','linewidth',2)
    ylim([0 5])
    if n==1
        xlabel('Position (cm)')
        ylabel('Firing rate (spike/sec)')
        legend('boxoff')
        set(gca,'ytick',[0 5])
    else
        set(gca,'xtick',[])
        set(gca,'ytick',[])
    end
    set(gca,'fontsize',16)
    box off
end

% Figure 1D - Example place cells in CA3:
figure('units','normalized','outerposition',[0.2 0.1 0.6 0.8])
for n=1:size(example_cells_spiking_activity_CA3,2)
    axes('position',[(n-1)*0.185+0.07 0.65 0.165 0.24])
    plot(bin_size*example_day_mouse_position_CA3,t_vec_CA3(end)-t_vec_CA3)
    chosen_cell_spikes=find(example_cells_spiking_activity_CA3(:,n)>0);
    chosen_cell_velocities=example_day_mouse_velocities_CA3(chosen_cell_spikes);
    
    num_spikes_chosen_cell=length(chosen_cell_spikes);
    for k=1:num_spikes_chosen_cell
        hold on
        if chosen_cell_velocities(k)>0
            plot(bin_size*example_day_mouse_position_CA3(chosen_cell_spikes(k)),t_vec_CA3(example_day_concatanation_index_CA3)-t_vec_CA3(chosen_cell_spikes(k)),'.','markersize',20,'color','r')
        else
            plot(bin_size*example_day_mouse_position_CA3(chosen_cell_spikes(k)),t_vec_CA3(example_day_concatanation_index_CA3)-t_vec_CA3(chosen_cell_spikes(k)),'.','markersize',20,'color',[0 0.8 0])
        end
        ylim([0 t_vec_CA3(example_day_concatanation_index_CA3)])
    end
    if n==1
        ylabel('Time (sec)')
    end
    set(gca,'ytick',[])
    set(gca,'xtick',[])
    hold on
    title(['Cell ' num2str(n)],'fontweight','normal')
    set(gca,'fontsize',16)
    if n==3
        text(50,800,'Figure 1D - CA3','fontsize',20,'fontweight','bold','HorizontalAlignment','Center')
    end
    
    axes('position',[(n-1)*0.185+0.07 0.39 0.165 0.24])
    plot(bin_size*example_day_mouse_position_CA3(example_day_concatanation_index_CA3+1:end),t_vec_CA3(example_day_concatanation_index_CA3)+t_vec_CA3(end)-t_vec_CA3(example_day_concatanation_index_CA3+1:end))
    chosen_cell_velocities=example_day_mouse_velocities_CA3(chosen_cell_spikes);
    
    num_spikes_chosen_cell=length(chosen_cell_spikes);
    for k=1:num_spikes_chosen_cell
        hold on
        if chosen_cell_velocities(k)>0
            plot(bin_size*example_day_mouse_position_CA3(chosen_cell_spikes(k)),t_vec_CA3(example_day_concatanation_index_CA3)+t_vec_CA3(end)-t_vec_CA3(chosen_cell_spikes(k)),'.','markersize',20,'color','r')
        else
            plot(bin_size*example_day_mouse_position_CA3(chosen_cell_spikes(k)),t_vec_CA3(example_day_concatanation_index_CA3)+t_vec_CA3(end)-t_vec_CA3(chosen_cell_spikes(k)),'.','markersize',20,'color',[0 0.8 0])
        end
        ylim([t_vec_CA3(example_day_concatanation_index_CA3+1) t_vec_CA3(end)])
    end
    if n==1
        ylabel('Time (sec)')
    end
    set(gca,'ytick',[])
    set(gca,'xtick',[])
    hold on
    set(gca,'fontsize',16)
    
    axes('position',[(n-1)*0.185+0.07 0.08 0.165 0.27])
    plot(-bin_size/2+bin_size*[1:total_num_bins],fs*(squeeze(example_day_tuning_curves_CA3(n,1,:)))/length(example_day_mouse_position_CA3)./squeeze(example_day_prior_distribution_CA3(:,1)),'color','r','linewidth',2)
    hold on
    plot(-bin_size/2+bin_size*[1:total_num_bins],fs*(squeeze(example_day_tuning_curves_CA3(n,2,:)))/length(example_day_mouse_position_CA3)./squeeze(example_day_prior_distribution_CA3(:,2)),'color','g','linewidth',2)
    ylim([0 5])
    if n==1
        xlabel('Position (cm)')
        ylabel('Firing rate (spike/sec)')
        legend('boxoff')
        set(gca,'ytick',[0 5])
    else
        set(gca,'xtick',[])
        set(gca,'ytick',[])
    end
    set(gca,'fontsize',16)
    box off
end


%% Figure 2 - CA3 place cells exhibit higher spatial tuning precision and short-term stability than CA1 place cells:

% Figure 2A-C - Spatial information:
SI_SSR_bit_spike_place_cells_novel_CA1=[];
SI_SSR_bit_spike_place_cells_novel_CA3=[];
for n=1:number_of_mice
    for k=1:maximal_number_of_sessions
        if mouse_group(n)==1
            if ~isempty(intersect(k,novel_sessions))
                SI_SSR_bit_spike_place_cells_novel_CA1=[SI_SSR_bit_spike_place_cells_novel_CA1,across_mice_data{n}.SI_SSR_bit_spike_sessions{k}];
            end
        elseif mouse_group(n)==2
            if ~isempty(intersect(k,novel_sessions))
                SI_SSR_bit_spike_place_cells_novel_CA3=[SI_SSR_bit_spike_place_cells_novel_CA3,across_mice_data{n}.SI_SSR_bit_spike_sessions{k}];
            end
        end
    end
end

SI_SSR_bit_spike_dynamics=nan(number_of_mice,maximal_number_of_sessions);
for n=1:number_of_mice
    SI_SSR_bit_spike_dynamics(n,:)=across_mice_data{n}.average_SI_SSR_bit_spike_sessions;
end

mean_SI_SSR_bit_spike_dynamics_CA1=mean(SI_SSR_bit_spike_dynamics(mouse_group==1,:),'omitnan');
std_SI_SSR_bit_spike_dynamics_CA1=std(SI_SSR_bit_spike_dynamics(mouse_group==1,:),'omitnan');
mean_SI_SSR_bit_spike_dynamics_CA3=mean(SI_SSR_bit_spike_dynamics(mouse_group==2,:),'omitnan');
std_SI_SSR_bit_spike_dynamics_CA3=std(SI_SSR_bit_spike_dynamics(mouse_group==2,:),'omitnan');
average_SI_SSR_bit_spike=mean(SI_SSR_bit_spike_dynamics(:,novel_sessions)');

% Figure 2A - Distribution of spatial information:
figure
x_vec=-0.5:0.2:4.9;
[n1,~]=hist(SI_SSR_bit_spike_place_cells_novel_CA1,x_vec);
[n2,~]=hist(SI_SSR_bit_spike_place_cells_novel_CA3,x_vec);
n1=n1./sum(n1);
n2=n2./sum(n2);
plot(x_vec,n1,'-b','linewidth',2);
hold on
plot(x_vec,n2,'-r','linewidth',2);
xlim([0 4])
ylim([0 0.2])
xlabel('Spatial information (bit/spike)')
ylabel('Fraction of cells')
legend('CA1','CA3')
legend boxoff
set(gca,'fontsize',16)
box off
axis square
title('Figure 2A')

% Figure 2B - Average spatial information:
figure
x_vec_1=1+0.3*(0.5-rand(1,size(average_SI_SSR_bit_spike(mouse_group==1),2)));
x_vec_2=2+0.3*(0.5-rand(1,size(average_SI_SSR_bit_spike(mouse_group==2),2)));
bar([mean(average_SI_SSR_bit_spike(mouse_group==1)),mean(average_SI_SSR_bit_spike(mouse_group==2))],0.5,'FaceColor','none')
hold on
errorbar([mean(average_SI_SSR_bit_spike(mouse_group==1)),mean(average_SI_SSR_bit_spike(mouse_group==2))],[std(average_SI_SSR_bit_spike(mouse_group==1))./sqrt(sum(mouse_group==1)),std(average_SI_SSR_bit_spike(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',3,'color','k')
hold on
scatter(x_vec_1,average_SI_SSR_bit_spike(mouse_group==1),15,'k','filled');
hold on
scatter(x_vec_2,average_SI_SSR_bit_spike(mouse_group==2),15,'k','filled');
xlim([0.5 2.5])
ylim([0 2])
set(gca,'ytick',0.5:0.5:2)
ylabel('Spatial information (bit/spike)')
set(gca,'fontsize',16)
box off
hold on
axis square
set(gca,'xtick',1:3)
set(gca,'xticklabels',{'CA1','CA3'})
title('Figure 2B','Fontweight','Bold')
set(gca,'fontsize',16)
set(gcf,'PaperPositionMode','auto')
axis square
box off

% Figure 2C - Spatial information over days:
figure
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_SI_SSR_bit_spike_dynamics_CA1(1:maximal_number_of_sessions/2),std_SI_SSR_bit_spike_dynamics_CA1(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
hold on
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_SI_SSR_bit_spike_dynamics_CA3(1:maximal_number_of_sessions/2),std_SI_SSR_bit_spike_dynamics_CA3(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
hold on
errorbar(days_vec(maximal_number_of_sessions/2+1:end),mean_SI_SSR_bit_spike_dynamics_CA1(maximal_number_of_sessions/2+1:end),std_SI_SSR_bit_spike_dynamics_CA1(maximal_number_of_sessions/2+1:end)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
hold on
errorbar(days_vec(maximal_number_of_sessions/2+1:end),mean_SI_SSR_bit_spike_dynamics_CA3(maximal_number_of_sessions/2+1:end),std_SI_SSR_bit_spike_dynamics_CA3(maximal_number_of_sessions/2+1:end)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
plot([17 17],[0 800],'--','color','k','linewidth',2)
set(gca,'xtick',1:10:31)
xlim([1 33])
ylim([0 2])
xlabel('Time (days)')
ylabel('Spatial information (bit/spike)')
title('Figure 2C','Fontweight','Bold')
set(gca,'fontsize',16)
legend('CA1','CA3','Location','Southwest')
legend('boxoff')
box off
axis square

% Figure 2D-I - PV correlations:
per_day_PV_correlations_bin_resolution=nan(number_of_mice,maximal_number_of_sessions,2*num_bins*2,2*num_bins*2);
for k=1:number_of_mice
    this_mouse_PV_correlations_bin_resolution=across_mice_data{k}.PV_correlations_bin_resolution;
    for n=1:maximal_number_of_sessions
        per_day_PV_correlations_bin_resolution(k,n,:,:)=this_mouse_PV_correlations_bin_resolution((n-1)*2*num_bins*2+1:n*2*num_bins*2,(n-1).*2*num_bins*2+1:n.*2*num_bins*2);
    end
end

per_day_PV_correlations_bin_resolution_CA1=squeeze(mean(per_day_PV_correlations_bin_resolution(mouse_group==1,:,:,:),'omitnan'));
per_day_PV_correlations_bin_resolution_CA3=squeeze(mean(per_day_PV_correlations_bin_resolution(mouse_group==2,:,:,:),'omitnan'));
PV_correlations_bin_resolution_novel_CA1=squeeze(mean(per_day_PV_correlations_bin_resolution_CA1(novel_sessions,:,:),'omitnan'));
PV_correlations_bin_resolution_novel_CA3=squeeze(mean(per_day_PV_correlations_bin_resolution_CA3(novel_sessions,:,:),'omitnan'));

within_day_PV_correlations=nan(number_of_mice,maximal_number_of_sessions);
for n=1:number_of_mice
    within_day_PV_correlations(n,:)=across_mice_data{n}.within_day_PV_correlations;
end
mean_within_day_PV_correlations_CA1=mean(within_day_PV_correlations(mouse_group==1,:),'omitnan');
std_within_day_PV_correlations_CA1=std(within_day_PV_correlations(mouse_group==1,:),'omitnan');
mean_within_day_PV_correlations_CA3=mean(within_day_PV_correlations(mouse_group==2,:),'omitnan');
std_within_day_PV_correlations_CA3=std(within_day_PV_correlations(mouse_group==2,:),'omitnan');

average_within_day_PV_correlations=mean(within_day_PV_correlations(:,novel_sessions)');

% PV correlation across directions:
PV_correlations_across_sides_maximal_shift=nan(number_of_mice,maximal_number_of_sessions);
possible_shifts=-5:5;
for n=1:number_of_mice
    for k=1:maximal_number_of_sessions
        this_day_PV_correlations=squeeze(per_day_PV_correlations_bin_resolution(n,k,:,:));
        
        temp_PV_correlations_vs_shift=nan(1,length(possible_shifts));
        for shift_ind=1:length(possible_shifts)
            temp_PV_correlations_this_shift=nan(1,2*num_bins);
            for s=1:num_bins-abs(possible_shifts(shift_ind))
                if possible_shifts(shift_ind)>=0
                    temp_PV_correlations_this_shift(s)=this_day_PV_correlations(s,s+num_bins+possible_shifts(shift_ind));
                    temp_PV_correlations_this_shift(num_bins+s)=this_day_PV_correlations(2*num_bins+s,2*num_bins+s+num_bins+possible_shifts(shift_ind));
                else
                    temp_PV_correlations_this_shift(s)=this_day_PV_correlations(s-possible_shifts(shift_ind),s+num_bins);
                    temp_PV_correlations_this_shift(num_bins+s)=this_day_PV_correlations(2*num_bins+s-possible_shifts(shift_ind),2*num_bins+s+num_bins);
                end
            end
            temp_PV_correlations_vs_shift(shift_ind)=mean(temp_PV_correlations_this_shift,'omitnan');
        end
        [this_maximal_shift_PV_correlation,~]=max(temp_PV_correlations_vs_shift);
        PV_correlations_across_sides_maximal_shift(n,k)=this_maximal_shift_PV_correlation;
    end
end

PV_correlations_novel_maximal_shift=mean(PV_correlations_across_sides_maximal_shift(:,novel_sessions),2,'omitnan');

mean_across_sides_maximal_shift_PV_dynamics_CA1=mean(PV_correlations_across_sides_maximal_shift(mouse_group==1,:),'omitnan');
std_across_sides_maximal_shift_PV_dynamics_CA1=std(PV_correlations_across_sides_maximal_shift(mouse_group==1,:),'omitnan');
mean_across_sides_maximal_shift_PV_dynamics_CA3=mean(PV_correlations_across_sides_maximal_shift(mouse_group==2,:),'omitnan');
std_across_sides_maximal_shift_PV_dynamics_CA3=std(PV_correlations_across_sides_maximal_shift(mouse_group==2,:),'omitnan');

PV_correlations_novel_maximal_shift_L=nan(1,number_of_mice);
PV_correlations_novel_maximal_shift_straight=nan(1,number_of_mice);
for n=1:number_of_mice
    if n>=3 && n<=5
        PV_correlations_novel_maximal_shift_L(n)=mean(PV_correlations_across_sides_maximal_shift(n,1:3),2,'omitnan');
        PV_correlations_novel_maximal_shift_straight(n)=mean(PV_correlations_across_sides_maximal_shift(n,9:11),2,'omitnan');
    else
        PV_correlations_novel_maximal_shift_L(n)=mean(PV_correlations_across_sides_maximal_shift(n,9:11),2,'omitnan');
        PV_correlations_novel_maximal_shift_straight(n)=mean(PV_correlations_across_sides_maximal_shift(n,1:3),2,'omitnan');
    end
end

% Figure 2D - Average PV correlations in CA1:
figure
imagesc(PV_correlations_bin_resolution_novel_CA1)
hold on
plot(0.5+[0 2*num_bins*2],[2*num_bins 2*num_bins]+0.5,'-w','linewidth',2)
plot([2*num_bins 2*num_bins]+0.5,0.5+[0 2*num_bins*2],'-w','linewidth',2)
plot(0.5+[0 2*num_bins*2],[num_bins num_bins]+0.5,'-k','linewidth',1)
plot(0.5+[0 2*num_bins*2],[3*num_bins 3*num_bins]+0.5,'-k','linewidth',1)
plot([num_bins num_bins]+0.5,0.5+[0 2*num_bins*2],'-k','linewidth',1)
plot([3*num_bins 3*num_bins]+0.5,0.5+[0 2*num_bins*2],'-k','linewidth',1)
cbh=colorbar;
set(cbh,'XTick',[-0.2 1])
colormap('jet')
axis square
set(gca,'xtick',[])
set(gca,'ytick',[])
box off
caxis([-0.2 1])
text(95,40,'PV correlation','fontsize',18,'fontweight','normal','rotation',90,'HorizontalAlignment','Center')
set(gca,'fontsize',16);
title('Figure 2D - CA1')

% Figure 2E - Average PV correlations in CA3:
figure
imagesc(PV_correlations_bin_resolution_novel_CA3)
hold on
plot(0.5+[0 2*num_bins*2],[2*num_bins 2*num_bins]+0.5,'-w','linewidth',2)
plot([2*num_bins 2*num_bins]+0.5,0.5+[0 2*num_bins*2],'-w','linewidth',2)
plot(0.5+[0 2*num_bins*2],[num_bins num_bins]+0.5,'-k','linewidth',1)
plot(0.5+[0 2*num_bins*2],[3*num_bins 3*num_bins]+0.5,'-k','linewidth',1)
plot([num_bins num_bins]+0.5,0.5+[0 2*num_bins*2],'-k','linewidth',1)
plot([3*num_bins 3*num_bins]+0.5,0.5+[0 2*num_bins*2],'-k','linewidth',1)
cbh=colorbar;
set(cbh,'XTick',[-0.2 1])
colormap('jet')
axis square
set(gca,'xtick',[])
set(gca,'ytick',[])
box off
caxis([-0.2 1])
text(95,40,'PV correlation','fontsize',18,'fontweight','normal','rotation',90,'HorizontalAlignment','Center')
set(gca,'fontsize',16);
title('Figure 2E - CA3')

% Figure 2F - Average PV correlation across directions:
figure
x_vec_1=1+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_2=2+0.3*(0.5-rand(1,sum(mouse_group==2)));
bar([1 2],[mean(PV_correlations_novel_maximal_shift(mouse_group==1)),mean(PV_correlations_novel_maximal_shift(mouse_group==2))],0.5,'FaceColor','none')
hold on
errorbar([1 2],[mean(PV_correlations_novel_maximal_shift(mouse_group==1)),mean(PV_correlations_novel_maximal_shift(mouse_group==2))],[std(PV_correlations_novel_maximal_shift(mouse_group==1))./sqrt(sum(mouse_group==1)),std(PV_correlations_novel_maximal_shift(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','k')
scatter(x_vec_1,PV_correlations_novel_maximal_shift(mouse_group==1),25,'k','filled');
scatter(x_vec_2,PV_correlations_novel_maximal_shift(mouse_group==2),25,'k','filled');
xlim([0.5 2.5])
ylim([-0.1 0.3])
ylabel('PV correlation acorss directions')
box off
hold on
axis square
set(gca,'xtick',[1 2])
set(gca,'xticklabels',{'CA1','CA3'})
set(gca,'fontsize',16)
axis square
box off
title ('Figure 2F - Within session')

% Figure 2G - PV correlation across directions over time:
figure
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_across_sides_maximal_shift_PV_dynamics_CA1(1:maximal_number_of_sessions/2),std_across_sides_maximal_shift_PV_dynamics_CA1(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
hold on
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_across_sides_maximal_shift_PV_dynamics_CA3(1:maximal_number_of_sessions/2),std_across_sides_maximal_shift_PV_dynamics_CA3(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
errorbar(days_vec(maximal_number_of_sessions/2+1:end),mean_across_sides_maximal_shift_PV_dynamics_CA1(maximal_number_of_sessions/2+1:end),std_across_sides_maximal_shift_PV_dynamics_CA1(maximal_number_of_sessions/2+1:end)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
errorbar(days_vec(maximal_number_of_sessions/2+1:end),mean_across_sides_maximal_shift_PV_dynamics_CA3(maximal_number_of_sessions/2+1:end),std_across_sides_maximal_shift_PV_dynamics_CA3(maximal_number_of_sessions/2+1:end)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
plot([17 17],[-0.1 1],'--','color','k','linewidth',2)
plot([1 33],[0 0],'--','color','k','linewidth',2)
set(gca,'xtick',1:10:31)
xlim([1 33])
ylim([-0.1 0.3])
set(gca,'ytick',-0.1:0.1:0.3)
xlabel('Time (days)')
ylabel('PV ocrrelations across directions')
set(gca,'fontsize',16)
legend('CA1','CA3')
legend boxoff
box off
axis square
title('Figure 2G - Within session')

% Figure 2G, inset - Comparison between environments of PV correlation across directions:
figure
x_vec_1=1+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_2=2+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_3=4+0.3*(0.5-rand(1,sum(mouse_group==2)));
x_vec_4=5+0.3*(0.5-rand(1,sum(mouse_group==2)));
plot([1 2],[mean(PV_correlations_novel_maximal_shift_L(mouse_group==1)), mean(PV_correlations_novel_maximal_shift_straight(mouse_group==1))],'.-b','markersize',30,'linewidth',2)
hold on
plot([0.5 5.5],[0 0],'--k','linewidth',2)
plot([4 5],[mean(PV_correlations_novel_maximal_shift_L(mouse_group==2)), mean(PV_correlations_novel_maximal_shift_straight(mouse_group==2))],'.-r','markersize',30,'linewidth',2)
errorbar([1 2],[mean(PV_correlations_novel_maximal_shift_L(mouse_group==1)),mean(PV_correlations_novel_maximal_shift_straight(mouse_group==1))],[std(PV_correlations_novel_maximal_shift_L(mouse_group==1))./sqrt(sum(mouse_group==1)),std(PV_correlations_novel_maximal_shift_straight(mouse_group==1))./sqrt(sum(mouse_group==1))],'.','linewidth',2,'color','b')
errorbar([4 5],[mean(PV_correlations_novel_maximal_shift_L(mouse_group==2)),mean(PV_correlations_novel_maximal_shift_straight(mouse_group==2))],[std(PV_correlations_novel_maximal_shift_L(mouse_group==2))./sqrt(sum(mouse_group==2)),std(PV_correlations_novel_maximal_shift_straight(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','r')
for n=1:sum(mouse_group==1)
    plot([x_vec_1(n), x_vec_2(n)],[PV_correlations_novel_maximal_shift_L(CA1_indexes(n)),PV_correlations_novel_maximal_shift_straight(CA1_indexes(n))],'-','color',[0.7 0.7 0.7])
end
for n=1:sum(mouse_group==2)
    plot([x_vec_3(n), x_vec_4(n)],[PV_correlations_novel_maximal_shift_L(CA3_indexes(n)),PV_correlations_novel_maximal_shift_straight(CA3_indexes(n))],'-','color',[0.7 0.7 0.7])
end
plot([1 2],[mean(PV_correlations_novel_maximal_shift_L(mouse_group==1)), mean(PV_correlations_novel_maximal_shift_straight(mouse_group==1))],'.-b','markersize',30,'linewidth',2)
plot([4 5],[mean(PV_correlations_novel_maximal_shift_L(mouse_group==2)), mean(PV_correlations_novel_maximal_shift_straight(mouse_group==2))],'.-r','markersize',30,'linewidth',2)
errorbar([1 2],[mean(PV_correlations_novel_maximal_shift_L(mouse_group==1)),mean(PV_correlations_novel_maximal_shift_straight(mouse_group==1))],[std(PV_correlations_novel_maximal_shift_L(mouse_group==1))./sqrt(sum(mouse_group==1)),std(PV_correlations_novel_maximal_shift_straight(mouse_group==1))./sqrt(sum(mouse_group==1))],'.','linewidth',2,'color','b')
errorbar([4 5],[mean(PV_correlations_novel_maximal_shift_L(mouse_group==2)),mean(PV_correlations_novel_maximal_shift_straight(mouse_group==2))],[std(PV_correlations_novel_maximal_shift_L(mouse_group==2))./sqrt(sum(mouse_group==2)),std(PV_correlations_novel_maximal_shift_straight(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','r')
xlim([0.5 5.5])
ylim([-0.1 0.3])
ylabel('PV correlation across directions')
box off
hold on
axis square
set(gca,'xtick',[1 2 4 5])
set(gca,'xticklabels',{'L','St' ,'L','St'})
set(gca,'fontsize',16)
axis square
box off
title('Figure 2G - inset')

% Figure 2H - Average within-day PV correlations:
figure
x_vec_1=1+0.3*(0.5-rand(1,size(average_within_day_PV_correlations(mouse_group==1),2)));
x_vec_2=2+0.3*(0.5-rand(1,size(average_within_day_PV_correlations(mouse_group==2),2)));
bar([mean(average_within_day_PV_correlations(mouse_group==1)),mean(average_within_day_PV_correlations(mouse_group==2))],0.5,'FaceColor','none')
hold on
errorbar([mean(average_within_day_PV_correlations(mouse_group==1)),mean(average_within_day_PV_correlations(mouse_group==2))],[std(average_within_day_PV_correlations(mouse_group==1))./sqrt(sum(mouse_group==1)),std(average_within_day_PV_correlations(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',3,'color','k')
scatter(x_vec_1,average_within_day_PV_correlations(mouse_group==1),15,'k','filled');
scatter(x_vec_2,average_within_day_PV_correlations(mouse_group==2),15,'k','filled');
xlim([0.5 2.5])
ylim([0 1])
set(gca,'ytick',0:0.2:1)
ylabel('PV correlation')
set(gca,'fontsize',16)
box off
hold on
axis square
set(gca,'xtick',1:3)
set(gca,'xticklabels',{'CA1','CA3'})
set(gca,'fontsize',16)
set(gcf,'PaperPositionMode','auto')
axis square
box off
title('Figure 2H - Within day')

% Figure 2I - Within-day PV correlations over time:
figure
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_within_day_PV_correlations_CA1(1:maximal_number_of_sessions/2),std_within_day_PV_correlations_CA1(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
hold on
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_within_day_PV_correlations_CA3(1:maximal_number_of_sessions/2),std_within_day_PV_correlations_CA3(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
errorbar(days_vec(maximal_number_of_sessions/2+1:end),mean_within_day_PV_correlations_CA1(maximal_number_of_sessions/2+1:end),std_within_day_PV_correlations_CA1(maximal_number_of_sessions/2+1:end)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
errorbar(days_vec(maximal_number_of_sessions/2+1:end),mean_within_day_PV_correlations_CA3(maximal_number_of_sessions/2+1:end),std_within_day_PV_correlations_CA3(maximal_number_of_sessions/2+1:end)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
plot([17 17],[0 800],'--','color','k','linewidth',2)
set(gca,'xtick',1:10:31)
xlim([1 33])
ylim([0 1])
set(gca,'ytick',-0.2:0.2:1)
xlabel('Time (days)')
ylabel('PV correlation')
set(gca,'fontsize',16)
box off
axis square
legend('CA1','CA3','Location','southeast')
legend('boxoff')
title('Figure 2I - Within day')

%% Figure 3 - CA3 exhibits more stable spatial representations over weeks than CA1:

% Figure 3A-D - Rate correlations and PV correlations across days:
PV_dynamics=nan(number_of_mice,maximal_number_of_sessions);
rate_dynamics=nan(number_of_mice,maximal_number_of_sessions);

PV_matrix=nan(number_of_mice,2*maximal_number_of_sessions,2*maximal_number_of_sessions);
rate_matrix=nan(number_of_mice,2*maximal_number_of_sessions,2*maximal_number_of_sessions);
maximal_PV_matrix=nan(number_of_mice,2*maximal_number_of_sessions,2*maximal_number_of_sessions);

for n=1:number_of_mice
    PV_dynamics(n,:)=across_mice_data{n}.mean_PV_dynamics(1:maximal_number_of_sessions);
    rate_dynamics(n,:)=across_mice_data{n}.mean_rate_dynamics(1:maximal_number_of_sessions);
    rate_matrix(n,:,:)=across_mice_data{n}.rate_correlations;
    PV_matrix(n,:,:)=across_mice_data{n}.PV_correlations;
    maximal_PV_matrix(n,:,:)=across_mice_data{n}.maximal_PV_correlations;
end

for n=1:number_of_mice
    for k=1:2*maximal_number_of_sessions
        PV_matrix(n,k,k)=1;
        maximal_PV_matrix(n,k,k)=1;
    end
end

mean_PV_dynamics_CA1=mean(PV_dynamics(mouse_group==1,:),'omitnan');
std_PV_dynamics_CA1=std(PV_dynamics(mouse_group==1,:),'omitnan');
mean_rate_dynamics_CA1=mean(rate_dynamics(mouse_group==1,:),'omitnan');
std_rate_dynamics_CA1=std(rate_dynamics(mouse_group==1,:),'omitnan');
average_PV_matrix_CA1=squeeze(mean(PV_matrix(mouse_group==1,:,:),'omitnan'));
average_rate_matrix_CA1=squeeze(mean(rate_matrix(mouse_group==1,:,:),'omitnan'));

mean_PV_dynamics_CA3=mean(PV_dynamics(mouse_group==2,:),'omitnan');
std_PV_dynamics_CA3=std(PV_dynamics(mouse_group==2,:),'omitnan');
mean_rate_dynamics_CA3=mean(rate_dynamics(mouse_group==2,:),'omitnan');
std_rate_dynamics_CA3=std(rate_dynamics(mouse_group==2,:),'omitnan');
average_PV_matrix_CA3=squeeze(mean(PV_matrix(mouse_group==2,:,:),'omitnan'));
average_rate_matrix_CA3=squeeze(mean(rate_matrix(mouse_group==2,:,:),'omitnan'));

averaged_across_environments_PV_matrix_CA1_temp=zeros(2,maximal_number_of_sessions*num_trials/2,maximal_number_of_sessions*num_trials/2);
averaged_across_environments_PV_matrix_CA1_temp(1,:,:)=average_PV_matrix_CA1(1:maximal_number_of_sessions*num_trials/2,1:maximal_number_of_sessions*num_trials/2);
averaged_across_environments_PV_matrix_CA1_temp(2,:,:)=average_PV_matrix_CA1(maximal_number_of_sessions*num_trials/2+1:end,maximal_number_of_sessions*num_trials/2+1:end);
averaged_across_environments_PV_matrix_CA1=squeeze(mean(averaged_across_environments_PV_matrix_CA1_temp,'omitnan'));
averaged_across_environments_PV_matrix_CA3_temp=zeros(2,maximal_number_of_sessions*num_trials/2,maximal_number_of_sessions*num_trials/2);
averaged_across_environments_PV_matrix_CA3_temp(1,:,:)=average_PV_matrix_CA3(1:maximal_number_of_sessions*num_trials/2,1:maximal_number_of_sessions*num_trials/2);
averaged_across_environments_PV_matrix_CA3_temp(2,:,:)=average_PV_matrix_CA3(maximal_number_of_sessions*num_trials/2+1:end,maximal_number_of_sessions*num_trials/2+1:end);
averaged_across_environments_PV_matrix_CA3=squeeze(mean(averaged_across_environments_PV_matrix_CA3_temp,'omitnan'));

averaged_across_environments_rate_matrix_CA1_temp=zeros(2,maximal_number_of_sessions*num_trials/2,maximal_number_of_sessions*num_trials/2);
averaged_across_environments_rate_matrix_CA1_temp(1,:,:)=average_rate_matrix_CA1(1:maximal_number_of_sessions*num_trials/2,1:maximal_number_of_sessions*num_trials/2);
averaged_across_environments_rate_matrix_CA1_temp(2,:,:)=average_rate_matrix_CA1(maximal_number_of_sessions*num_trials/2+1:end,maximal_number_of_sessions*num_trials/2+1:end);
averaged_across_environments_rate_matrix_CA1=squeeze(mean(averaged_across_environments_rate_matrix_CA1_temp,'omitnan'));
averaged_across_environments_rate_matrix_CA3_temp=zeros(2,maximal_number_of_sessions*num_trials/2,maximal_number_of_sessions*num_trials/2);
averaged_across_environments_rate_matrix_CA3_temp(1,:,:)=average_rate_matrix_CA3(1:maximal_number_of_sessions*num_trials/2,1:maximal_number_of_sessions*num_trials/2);
averaged_across_environments_rate_matrix_CA3_temp(2,:,:)=average_rate_matrix_CA3(maximal_number_of_sessions*num_trials/2+1:end,maximal_number_of_sessions*num_trials/2+1:end);
averaged_across_environments_rate_matrix_CA3=squeeze(mean(averaged_across_environments_rate_matrix_CA3_temp,'omitnan'));

% Figure 3A - Rate correlation matrix across days:
figure('units','normalized','outerposition',[0.2 0.2 0.6 0.6])
axes('position',[0.06 0.15 0.35 0.75])
imagesc(averaged_across_environments_rate_matrix_CA1)
title('CA1')
xlabel('Day')
ylabel('Day')
caxis([0 1])
set(gca,'xtick',0.5+ceil(num_trials/2):num_trials:maximal_number_of_sessions*num_trials)
set(gca,'xticklabel',days_vec(1:end))
set(gca,'ytick',0.5+ceil(num_trials/2):num_trials:maximal_number_of_sessions*num_trials)
set(gca,'yticklabel',days_vec(1:end))
set(gca,'fontsize',16)
axis square
colormap('jet')
for n=1:maximal_number_of_sessions/2-1
    hold on
    plot(0.5+[0 maximal_number_of_sessions*num_trials],0.5+[n*num_trials n*num_trials],'linewidth',2,'color','k')
    hold on
    plot(0.5+[n*num_trials n*num_trials],0.5+[0 maximal_number_of_sessions*num_trials],'linewidth',2,'color','k')
end
axes('position',[0.42 0.2 0.015 0.65])
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim([0 1])
ylim([0 1])
color_map=colormap('jet');
for n=1:size(color_map,1)
    hold on
    p=patch([0 1 1 0],[n/size(color_map,1) n/size(color_map,1) (n-1)/size(color_map,1) (n-1)/size(color_map,1)],color_map(n,:));
    set(p,'FaceAlpha',1,'EdgeColor','none');
end
plot([0 1 1 0 0],[0 0 1 1 0],'color','k')
text(2,0.5,'Rate correlation','fontsize',16,'fontweight','normal','rotation',90,'HorizontalAlignment','Center')
text(1.2,0,'0','fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
text(1.2,1,'1','fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
set(gca,'fontsize',16)
text(1.5,1.15,'Figure 3A','fontsize',20,'fontweight','bold','HorizontalAlignment','Left')

axes('position',[0.58 0.15 0.35 0.75])
imagesc(averaged_across_environments_rate_matrix_CA3)
title('CA3')
xlabel('Day')
ylabel('Day')
caxis([0 1])
set(gca,'xtick',0.5+ceil(num_trials/2):num_trials:maximal_number_of_sessions*num_trials)
set(gca,'xticklabel',days_vec(1:end))
set(gca,'ytick',0.5+ceil(num_trials/2):num_trials:maximal_number_of_sessions*num_trials)
set(gca,'yticklabel',days_vec(1:end))
set(gca,'fontsize',16)
axis square
colormap('jet')
for n=1:maximal_number_of_sessions/2-1
    hold on
    plot(0.5+[0 maximal_number_of_sessions*num_trials],0.5+[n*num_trials n*num_trials],'linewidth',2,'color','k')
    hold on
    plot(0.5+[n*num_trials n*num_trials],0.5+[0 maximal_number_of_sessions*num_trials],'linewidth',2,'color','k')
end
axes('position',[0.94 0.2 0.015 0.65])
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim([0 1])
ylim([0 1])
color_map=colormap('jet');
for n=1:size(color_map,1)
    hold on
    p=patch([0 1 1 0],[n/size(color_map,1) n/size(color_map,1) (n-1)/size(color_map,1) (n-1)/size(color_map,1)],color_map(n,:));
    set(p,'FaceAlpha',1,'EdgeColor','none');
end
plot([0 1 1 0 0],[0 0 1 1 0],'color','k')
text(2,0.5,'Rate correlation','fontsize',16,'fontweight','normal','rotation',90,'HorizontalAlignment','Center')
text(1.2,0,'0','fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
text(1.2,1,'1','fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
set(gca,'fontsize',16)

% Figure 3B - Rate correlation versus elapsed time:
figure
errorbar(elapsed_days_vec(2:end),mean_rate_dynamics_CA1(2:end),std_rate_dynamics_CA1(2:end)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
hold on
errorbar(elapsed_days_vec(2:end),mean_rate_dynamics_CA3(2:end),std_rate_dynamics_CA3(2:end)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
errorbar(elapsed_days_vec(1),mean_rate_dynamics_CA1(1),std_rate_dynamics_CA1(1)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
errorbar(elapsed_days_vec(1),mean_rate_dynamics_CA3(1),std_rate_dynamics_CA3(1)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
xlim([-0.5 14])
ylim([0 1])
xlabel('Elapsed time (days)')
ylabel('Rate correlation')
set(gca,'fontsize',16)
legend('CA1','CA3')
legend('boxoff')
box off
axis square
title('Figure 3B - Across days')

% Figure 3C - PV correlation matrix across days:
figure('units','normalized','outerposition',[0.2 0.2 0.6 0.6])
axes('position',[0.06 0.15 0.35 0.75])
imagesc(averaged_across_environments_PV_matrix_CA1)
title('CA1')
xlabel('Day')
ylabel('Day')
caxis([0 1])
set(gca,'xtick',0.5+ceil(num_trials/2):num_trials:maximal_number_of_sessions*num_trials)
set(gca,'xticklabel',days_vec(1:end))
set(gca,'ytick',0.5+ceil(num_trials/2):num_trials:maximal_number_of_sessions*num_trials)
set(gca,'yticklabel',days_vec(1:end))
set(gca,'fontsize',16)
axis square
colormap('jet')
for n=1:maximal_number_of_sessions/2-1
    hold on
    plot(0.5+[0 maximal_number_of_sessions*num_trials],0.5+[n*num_trials n*num_trials],'linewidth',2,'color','k')
    hold on
    plot(0.5+[n*num_trials n*num_trials],0.5+[0 maximal_number_of_sessions*num_trials],'linewidth',2,'color','k')
end
axes('position',[0.42 0.2 0.015 0.65])
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim([0 1])
ylim([0 1])
color_map=colormap('jet');
for n=1:size(color_map,1)
    hold on
    p=patch([0 1 1 0],[n/size(color_map,1) n/size(color_map,1) (n-1)/size(color_map,1) (n-1)/size(color_map,1)],color_map(n,:));
    set(p,'FaceAlpha',1,'EdgeColor','none');
end
plot([0 1 1 0 0],[0 0 1 1 0],'color','k')
text(2,0.5,'PV correlation','fontsize',16,'fontweight','normal','rotation',90,'HorizontalAlignment','Center')
text(1.2,0,'0','fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
text(1.2,1,'1','fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
set(gca,'fontsize',16)
text(1.5,1.15,'Figure 3C','fontsize',20,'fontweight','bold','HorizontalAlignment','Left')

axes('position',[0.58 0.15 0.35 0.75])
imagesc(averaged_across_environments_PV_matrix_CA3)
title('CA3')
xlabel('Day')
ylabel('Day')
caxis([0 1])
set(gca,'xtick',0.5+ceil(num_trials/2):num_trials:maximal_number_of_sessions*num_trials)
set(gca,'xticklabel',days_vec(1:end))
set(gca,'ytick',0.5+ceil(num_trials/2):num_trials:maximal_number_of_sessions*num_trials)
set(gca,'yticklabel',days_vec(1:end))
set(gca,'fontsize',16)
axis square
colormap('jet')
for n=1:maximal_number_of_sessions/2-1
    hold on
    plot(0.5+[0 maximal_number_of_sessions*num_trials],0.5+[n*num_trials n*num_trials],'linewidth',2,'color','k')
    hold on
    plot(0.5+[n*num_trials n*num_trials],0.5+[0 maximal_number_of_sessions*num_trials],'linewidth',2,'color','k')
end
axes('position',[0.94 0.2 0.015 0.65])
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim([0 1])
ylim([0 1])
color_map=colormap('jet');
for n=1:size(color_map,1)
    hold on
    p=patch([0 1 1 0],[n/size(color_map,1) n/size(color_map,1) (n-1)/size(color_map,1) (n-1)/size(color_map,1)],color_map(n,:));
    set(p,'FaceAlpha',1,'EdgeColor','none');
end
plot([0 1 1 0 0],[0 0 1 1 0],'color','k')
text(2,0.5,'PV correlation','fontsize',16,'fontweight','normal','rotation',90,'HorizontalAlignment','Center')
text(1.2,0,'0','fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
text(1.2,1,'1','fontsize',16,'fontweight','normal','HorizontalAlignment','Left')
set(gca,'fontsize',16)

% Figure 3D - PV correlation versus elapsed time:
figure
errorbar(elapsed_days_vec(2:end),mean_PV_dynamics_CA1(2:end),std_PV_dynamics_CA1(2:end)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
hold on
errorbar(elapsed_days_vec(2:end),mean_PV_dynamics_CA3(2:end),std_PV_dynamics_CA3(2:end)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
errorbar(elapsed_days_vec(1),mean_PV_dynamics_CA1(1),std_PV_dynamics_CA1(1)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
errorbar(elapsed_days_vec(1),mean_PV_dynamics_CA3(1),std_PV_dynamics_CA3(1)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
xlim([-0.5 14])
ylim([0 1])
set(gca,'ytick',-0.2:0.2:1)
xlabel('Elapsed time (days)')
ylabel('PV correlation')
set(gca,'fontsize',16)
box off
axis square
legend('CA1','CA3')
legend('boxoff')
title('Figure 3D - Across days')

%% Figure 4 - Hippocampal CA3 is organized into functionally related place cell assemblies:

% Loading all the necessary data for Figure 4:
all_within_day_pairwise_correlations_per_mouse_novel=cell(1,number_of_mice);
same_bin_within_day_pairwise_correlations_per_mouse_novel=cell(1,number_of_mice);
same_bin_across_sides_pairwise_correlations_per_mouse_novel=cell(1,number_of_mice);
diff_bins_within_day_pairwise_correlations_per_mouse_novel=cell(1,number_of_mice);
diff_bins_across_sides_pairwise_correlations_per_mouse_novel=cell(1,number_of_mice);
for n=1:number_of_mice
    all_within_day_pairwise_correlations_per_mouse_novel{n}=across_mice_data{n}.all_within_day_pairwise_correlations_per_mouse_novel;
    same_bin_within_day_pairwise_correlations_per_mouse_novel{n}=across_mice_data{n}.same_bin_within_day_pairwise_correlations_per_mouse_novel;
    same_bin_across_sides_pairwise_correlations_per_mouse_novel{n}=across_mice_data{n}.same_bin_across_sides_pairwise_correlations_per_mouse_novel;
    diff_bins_within_day_pairwise_correlations_per_mouse_novel{n}=across_mice_data{n}.diff_bins_within_day_pairwise_correlations_per_mouse_novel;
    diff_bins_across_sides_pairwise_correlations_per_mouse_novel{n}=across_mice_data{n}.diff_bins_across_sides_pairwise_correlations_per_mouse_novel;
end

% Distribution of all pairwise correlations:
corr_vec=-0.965:0.03:0.985;
all_within_day_pairwise_correlations_CA1_novel=[];
all_within_day_pairwise_correlations_CA3_novel=[];
for n=1:number_of_mice
    if mouse_group(n)==1
        all_within_day_pairwise_correlations_CA1_novel=[all_within_day_pairwise_correlations_CA1_novel,all_within_day_pairwise_correlations_per_mouse_novel{n}];
    elseif mouse_group(n)==2
        all_within_day_pairwise_correlations_CA3_novel=[all_within_day_pairwise_correlations_CA3_novel,all_within_day_pairwise_correlations_per_mouse_novel{n}];
    end
end
within_day_dist_pairwise_CA1_novel=hist(all_within_day_pairwise_correlations_CA1_novel,corr_vec);
normalized_within_day_dist_pairwise_CA1_novel=within_day_dist_pairwise_CA1_novel./sum(within_day_dist_pairwise_CA1_novel);
within_day_dist_pairwise_CA3_novel=hist(all_within_day_pairwise_correlations_CA3_novel,corr_vec);
normalized_within_day_dist_pairwise_CA3_novel=within_day_dist_pairwise_CA3_novel./sum(within_day_dist_pairwise_CA3_novel);

high_corr_threshold=0.7;
high_corr_index=find(corr_vec>high_corr_threshold);
normalized_all_within_day_pair_corr_per_mouse_novel=nan(number_of_mice,length(corr_vec));
for n=1:number_of_mice
    temp_hist=hist(all_within_day_pairwise_correlations_per_mouse_novel{n},corr_vec);
    normalized_temp_hist=temp_hist./sum(temp_hist);
    normalized_all_within_day_pair_corr_per_mouse_novel(n,:)=normalized_temp_hist;
end
fraction_pairs_high_corr_all_cells=sum(normalized_all_within_day_pair_corr_per_mouse_novel(:,high_corr_index)');
fraction_pairs_high_corr_all_cells_CA1=fraction_pairs_high_corr_all_cells(mouse_group==1);
fraction_pairs_high_corr_all_cells_CA3=fraction_pairs_high_corr_all_cells(mouse_group==2);

% Figure 4B - Distribution of all pairwise correlations:
figure
plot(corr_vec,normalized_within_day_dist_pairwise_CA1_novel,'-b','linewidth',2);
hold on
plot(corr_vec,normalized_within_day_dist_pairwise_CA3_novel,'-r','linewidth',2);
xlim([-0.5 1])
ylim([0 0.12])
xlabel('Tuning curve correlation')
ylabel('Fraction of cell-pairs')
set(gca,'xtick',-0.5:0.5:1)
legend('CA1','CA3')
legend boxoff
axis square
set(gca,'fontsize',16)
box off
title('Figure 4B - All positions')

% Figure 4B, inset - Fractoin of pairwise correlations > 0.7:
figure
x_vec_1=1+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_2=2+0.3*(0.5-rand(1,sum(mouse_group==2)));
bar([1 2],[mean(fraction_pairs_high_corr_all_cells_CA1),mean(fraction_pairs_high_corr_all_cells_CA3)],0.5,'FaceColor','none')
hold on
errorbar([1 2],[mean(fraction_pairs_high_corr_all_cells_CA1),mean(fraction_pairs_high_corr_all_cells_CA3)],[std(fraction_pairs_high_corr_all_cells_CA1)./sqrt(sum(mouse_group==1)),std(fraction_pairs_high_corr_all_cells_CA3)./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','k')
scatter(x_vec_1,fraction_pairs_high_corr_all_cells_CA1,25,'k','filled');
scatter(x_vec_2,fraction_pairs_high_corr_all_cells_CA3,25,'k','filled');
xlim([0.5 2.5])
ylim([0 0.1])
ylabel('Fraction of cell pairs')
box off
hold on
axis square
set(gca,'xtick',[1 2])
set(gca,'xticklabels',{'CA1','CA3'})
set(gca,'fontsize',16)
axis square
box off
title('Figure 4B - inset')

% Distribution of pairwise correlations between place cells with the same preferred position:
corr_vec=-0.99:0.02:0.99;
same_bin_within_day_pairwise_correlations_CA1_novel=[];
same_bin_within_day_pairwise_correlations_CA3_novel=[];
for n=1:number_of_mice
    if mouse_group(n)==1
        same_bin_within_day_pairwise_correlations_CA1_novel=[same_bin_within_day_pairwise_correlations_CA1_novel,same_bin_within_day_pairwise_correlations_per_mouse_novel{n}];
    elseif mouse_group(n)==2
        same_bin_within_day_pairwise_correlations_CA3_novel=[same_bin_within_day_pairwise_correlations_CA3_novel,same_bin_within_day_pairwise_correlations_per_mouse_novel{n}];
    end
end
same_bin_within_day_pairwise_dist_CA1_novel=hist(same_bin_within_day_pairwise_correlations_CA1_novel,corr_vec);
normalized_same_bin_within_day_pairwise_dist_CA1_novel=same_bin_within_day_pairwise_dist_CA1_novel./sum(same_bin_within_day_pairwise_dist_CA1_novel);
same_bin_within_day_pairwise_dist_CA3_novel=hist(same_bin_within_day_pairwise_correlations_CA3_novel,corr_vec);
normalized_same_bin_within_day_pairwise_dist_CA3_novel=same_bin_within_day_pairwise_dist_CA3_novel./sum(same_bin_within_day_pairwise_dist_CA3_novel);

normalized_same_bin_within_day_pair_corr_per_mouse_novel=nan(number_of_mice,length(corr_vec));
for n=1:number_of_mice
    temp_hist=hist(same_bin_within_day_pairwise_correlations_per_mouse_novel{n},corr_vec);
    normalized_temp_hist=temp_hist./sum(temp_hist);
    normalized_same_bin_within_day_pair_corr_per_mouse_novel(n,:)=normalized_temp_hist;
end

high_corr_threshold=0.7;
high_corr_index=find(corr_vec>high_corr_threshold);
fraction_pairs_high_corr_all_cells=sum(normalized_same_bin_within_day_pair_corr_per_mouse_novel(:,high_corr_index)');
fraction_pairs_high_corr_all_cells_CA1=fraction_pairs_high_corr_all_cells(mouse_group==1);
fraction_pairs_high_corr_all_cells_CA3=fraction_pairs_high_corr_all_cells(mouse_group==2);

% Figure 4D - Distribution of pairwise correlations between place cells with the same preferred position:
figure
plot(corr_vec,normalized_same_bin_within_day_pairwise_dist_CA1_novel,'-b','linewidth',2)
hold on
plot(corr_vec,normalized_same_bin_within_day_pairwise_dist_CA3_novel,'-r','linewidth',2)
xlim([-0.5 1])
ylim([0 0.05])
xlabel('Tuning curve correlation')
ylabel('Fraction of cell-pairs')
set(gca,'xtick',-0.5:0.5:1)
legend('CA1','CA3','Location','Northwest')
legend boxoff
axis square
set(gca,'fontsize',16)
box off
title('Figure 4D - Same position')

% Figure 4D, inset - Fraction of pairwise correlations > 0.7 between place cells with the same preferred position:
figure
x_vec_1=1+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_2=2+0.3*(0.5-rand(1,sum(mouse_group==2)));
bar([1 2],[mean(fraction_pairs_high_corr_all_cells_CA1),mean(fraction_pairs_high_corr_all_cells_CA3)],0.5,'FaceColor','none')
hold on
errorbar([1 2],[mean(fraction_pairs_high_corr_all_cells_CA1),mean(fraction_pairs_high_corr_all_cells_CA3)],[std(fraction_pairs_high_corr_all_cells_CA1)./sqrt(sum(mouse_group==1)),std(fraction_pairs_high_corr_all_cells_CA3)./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','k')
scatter(x_vec_1,fraction_pairs_high_corr_all_cells_CA1,25,'k','filled');
scatter(x_vec_2,fraction_pairs_high_corr_all_cells_CA3,25,'k','filled');
xlim([0.5 2.5])
ylim([0 1])
ylabel('Fraction of cell pairs')
box off
hold on
axis square
set(gca,'xtick',[1 2])
set(gca,'xticklabels',{'CA1','CA3'})
set(gca,'fontsize',16)
axis square
box off
title('Figure 4D - inset')

% Figure 4F-G - Comparing within-direction against across-directions
% distribution of pairwise correlations between place cells with the same
% preferred position:
corr_vec=-0.975:0.05:0.975;
normalized_same_bin_within_day_pair_corr_per_mouse_novel=nan(number_of_mice,length(corr_vec));
normalized_same_bin_across_sides_pair_corr_per_mouse_novel=nan(number_of_mice,length(corr_vec));
for n=1:number_of_mice
    temp_hist=hist(same_bin_within_day_pairwise_correlations_per_mouse_novel{n},corr_vec);
    normalized_temp_hist=temp_hist./sum(temp_hist);
    normalized_same_bin_within_day_pair_corr_per_mouse_novel(n,:)=normalized_temp_hist;
    
    temp_hist=hist(same_bin_across_sides_pairwise_correlations_per_mouse_novel{n},corr_vec);
    normalized_temp_hist=temp_hist./sum(temp_hist);
    normalized_same_bin_across_sides_pair_corr_per_mouse_novel(n,:)=normalized_temp_hist;
end

cdf_same_bin_within_day_pair_corr_per_mouse_novel=cumsum(normalized_same_bin_within_day_pair_corr_per_mouse_novel')';
cdf_same_bin_across_sides_pair_corr_per_mouse_novel=cumsum(normalized_same_bin_across_sides_pair_corr_per_mouse_novel')';

cdf_diff_same_bin_within_day_pair_corr_CA1_novel=(1-fliplr(cdf_same_bin_within_day_pair_corr_per_mouse_novel(mouse_group==1,:)))-(1-fliplr(cdf_same_bin_across_sides_pair_corr_per_mouse_novel(mouse_group==1,:)));
mean_cdf_diff_same_bin_within_day_pair_corr_CA1_novel=mean(cdf_diff_same_bin_within_day_pair_corr_CA1_novel);
sd_cdf_diff_same_bin_within_day_pair_corr_CA1_novel=std(cdf_diff_same_bin_within_day_pair_corr_CA1_novel);

cdf_diff_same_bin_within_day_pair_corr_CA3_novel=(1-fliplr(cdf_same_bin_within_day_pair_corr_per_mouse_novel(mouse_group==2,:)))-(1-fliplr(cdf_same_bin_across_sides_pair_corr_per_mouse_novel(mouse_group==2,:)));
mean_cdf_diff_same_bin_within_day_pair_corr_CA3_novel=mean(cdf_diff_same_bin_within_day_pair_corr_CA3_novel,'omitnan');
sd_cdf_diff_same_bin_within_day_pair_corr_CA3_novel=std(cdf_diff_same_bin_within_day_pair_corr_CA3_novel,'omitnan');

[~,per_mouse_max_diff_index_CA1_novel]=max(abs(cdf_diff_same_bin_within_day_pair_corr_CA1_novel'));
[~,per_mouse_max_diff_index_CA3_novel]=max(abs(cdf_diff_same_bin_within_day_pair_corr_CA3_novel'));
max_diff_CA1_novel=nan(1,sum(mouse_group==1));
for n=1:sum(mouse_group==1)
    max_diff_CA1_novel(n)=cdf_diff_same_bin_within_day_pair_corr_CA1_novel(n,per_mouse_max_diff_index_CA1_novel(n));
end
max_diff_CA3_novel=nan(1,sum(mouse_group==2));
for n=1:sum(mouse_group==2)
    max_diff_CA3_novel(n)=cdf_diff_same_bin_within_day_pair_corr_CA3_novel(n,per_mouse_max_diff_index_CA3_novel(n));
end

% Figure 4F - Cumulative distributions of within-direction and across-directions pairwise correlations:
figure
errorbar(corr_vec,1-fliplr(mean(cdf_same_bin_within_day_pair_corr_per_mouse_novel(mouse_group==1,:))),fliplr(std(cdf_same_bin_within_day_pair_corr_per_mouse_novel(mouse_group==1,:)))./sqrt(sum(mouse_group==1)),'-ob','markersize',5,'linewidth',1)
hold on
plot(corr_vec,1-fliplr(mean(cdf_same_bin_across_sides_pair_corr_per_mouse_novel(mouse_group==1,:))),'--b','linewidth',2)
errorbar(corr_vec,1-fliplr(mean(cdf_same_bin_within_day_pair_corr_per_mouse_novel(mouse_group==2,:))),fliplr(std(cdf_same_bin_within_day_pair_corr_per_mouse_novel(mouse_group==2,:)))./sqrt(sum(mouse_group==2)),'-or','markersize',5,'linewidth',1)
plot(corr_vec,1-fliplr(mean(cdf_same_bin_across_sides_pair_corr_per_mouse_novel(mouse_group==2,:))),'--r','linewidth',2)
xlim([-1 0])
ylim([0 1])
xlabel('Tuning curve correlation')
ylabel('Cumulative fraction of cell-pairs')
set(gca,'xtick',-1:0.5:1)
set(gca,'xticklabels',1:-0.5:-1)
legend('CA1 - Same direction','CA1 - Across directions','CA3 - Same direction','CA3 - Across directions','Location','Southeast')
legend boxoff
axis square
set(gca,'fontsize',16)
box off
title('Figure 4F - Same position')

% Figure 4F - Difference in cumulative distributions between the within-direction and across-directions pairwise correlations:
figure
errorbar(corr_vec,mean_cdf_diff_same_bin_within_day_pair_corr_CA1_novel,sd_cdf_diff_same_bin_within_day_pair_corr_CA1_novel./sqrt(sum(mouse_group==1)),'-b','linewidth',2)
hold on
errorbar(corr_vec,mean_cdf_diff_same_bin_within_day_pair_corr_CA3_novel,sd_cdf_diff_same_bin_within_day_pair_corr_CA3_novel./sqrt(sum(mouse_group==2)),'-r','linewidth',2)
plot([-1 0],[0 0],'--k','linewidth',2)
xlim([-1 0])
ylim([-0.05 0.1])
xlabel('Tuning curve correlation')
ylabel('Cumulative fraction difference')
set(gca,'xtick',-1:0.5:1)
set(gca,'xticklabels',1:-0.5:-1)
legend('CA1','CA3')
legend boxoff
axis square
set(gca,'fontsize',16)
box off
title('Novel sessions')
title('Figure 4G - Same position')

% Figure 4G, inset - Maximal difference in cumulative distributions between the within-direction and across-directions pairwise correlations:
figure
x_vec_1=1+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_2=2+0.3*(0.5-rand(1,sum(mouse_group==2)));
bar([1 2],[mean(max_diff_CA1_novel),mean(max_diff_CA3_novel)],0.5,'FaceColor','none')
hold on
errorbar([1 2],[mean(max_diff_CA1_novel),mean(max_diff_CA3_novel)],[std(max_diff_CA1_novel)./sqrt(sum(mouse_group==1)),std(max_diff_CA3_novel)./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','k')
scatter(x_vec_1,max_diff_CA1_novel,25,'k','filled');
scatter(x_vec_2,max_diff_CA3_novel,25,'k','filled');
xlim([0.5 2.5])
ylim([-0.05 0.1])
ylabel('Cumulative fraction difference')
box off
hold on
axis square
set(gca,'xtick',[1 2])
set(gca,'xticklabels',{'CA1','CA3'})
set(gca,'fontsize',16)
axis square
box off
title('Figure 4G - inset')


% Figure 4H-I - Comparing within-direction against across-directions
% distribution of pairwise correlations between place cells with the
% different preferred positions:
corr_vec=-0.975:0.05:0.975;
diff_bins_vec=0:1:5;
normalized_diff_bins_within_day_pair_corr_per_mouse_novel=nan(number_of_mice,length(diff_bins_vec),length(corr_vec));
for n=1:number_of_mice
    for k=1:length(diff_bins_vec)
        temp_hist=hist(diff_bins_within_day_pairwise_correlations_per_mouse_novel{n}{k},corr_vec);
        normalized_temp_hist=temp_hist./sum(temp_hist);
        normalized_diff_bins_within_day_pair_corr_per_mouse_novel(n,k,:)=normalized_temp_hist;
    end
end

normalized_diff_bins_across_sides_pair_corr_per_mouse_novel=nan(number_of_mice,length(diff_bins_vec),length(corr_vec));
for n=1:number_of_mice
    for k=1:length(diff_bins_vec)
        temp_hist=hist(diff_bins_across_sides_pairwise_correlations_per_mouse_novel{n}{k},corr_vec);
        normalized_temp_hist=temp_hist./sum(temp_hist);
        normalized_diff_bins_across_sides_pair_corr_per_mouse_novel(n,k,:)=normalized_temp_hist;
    end
end

cdf_diff_bins_within_day_pair_corr_per_mouse_novel=nan(number_of_mice,length(diff_bins_vec),length(corr_vec));
for n=1:number_of_mice
    for k=1:length(diff_bins_vec)
        cdf_diff_bins_within_day_pair_corr_per_mouse_novel(n,k,:)=cumsum(squeeze(normalized_diff_bins_within_day_pair_corr_per_mouse_novel(n,k,:))')';
    end
end

cdf_diff_bins_across_sides_pair_corr_per_mouse_novel=nan(number_of_mice,length(diff_bins_vec),length(corr_vec));
for n=1:number_of_mice
    for k=1:length(diff_bins_vec)
        cdf_diff_bins_across_sides_pair_corr_per_mouse_novel(n,k,:)=cumsum(squeeze(normalized_diff_bins_across_sides_pair_corr_per_mouse_novel(n,k,:))')';
    end
end

cdf_diff_diff_bins_within_day_pair_corr_CA1_novel=nan(sum(mouse_group==1),length(diff_bins_vec),length(corr_vec));
for n=1:sum(mouse_group==1)
    for k=1:length(diff_bins_vec)
        cdf_diff_diff_bins_within_day_pair_corr_CA1_novel(n,k,:)=(1-fliplr(squeeze(cdf_diff_bins_within_day_pair_corr_per_mouse_novel(CA1_indexes(n),k,:))))-(1-fliplr(squeeze(cdf_diff_bins_across_sides_pair_corr_per_mouse_novel(CA1_indexes(n),k,:))));
    end
end

cdf_diff_diff_bins_within_day_pair_corr_CA3_novel=nan(sum(mouse_group==1),length(diff_bins_vec),length(corr_vec));
for n=1:sum(mouse_group==2)
    for k=1:length(diff_bins_vec)
        cdf_diff_diff_bins_within_day_pair_corr_CA3_novel(n,k,:)=(1-fliplr(squeeze(cdf_diff_bins_within_day_pair_corr_per_mouse_novel(CA3_indexes(n),k,:))))-(1-fliplr(squeeze(cdf_diff_bins_across_sides_pair_corr_per_mouse_novel(CA3_indexes(n),k,:))));
    end
end

mean_cdf_diff_diff_bins_within_day_pair_corr_CA1_novel=squeeze(mean(cdf_diff_diff_bins_within_day_pair_corr_CA1_novel));
sd_cdf_diff_diff_bins_within_day_pair_corr_CA1_novel=squeeze(std(cdf_diff_diff_bins_within_day_pair_corr_CA1_novel));
mean_cdf_diff_diff_bins_within_day_pair_corr_CA3_novel=squeeze(mean(cdf_diff_diff_bins_within_day_pair_corr_CA3_novel));
sd_cdf_diff_diff_bins_within_day_pair_corr_CA3_novel=squeeze(std(cdf_diff_diff_bins_within_day_pair_corr_CA3_novel));

% Figure 4H - Difference in cumulative distributions for cells with different preferred positions (CA1):
smoothing_kernel=gausswin(5);
figure
color_vec_CA1=[0 1 1 ; 0 0.66 1; 0 0.33 1; 0 0 0.8 ; 0 0 0.5; 0 0 0];
hold on
for n=1:6
    plot(corr_vec,conv(fliplr(mean_cdf_diff_diff_bins_within_day_pair_corr_CA1_novel(n,:)),smoothing_kernel,'same')./conv(ones(1,length(corr_vec)),smoothing_kernel,'same'),'-','linewidth',2,'color',color_vec_CA1(n,:))
end
plot([-1 0],[0 0],'--k','linewidth',2)
xlim([-1 0.5])
ylim([-0.05 0.1])
xlabel('Tuning curve correlation')
ylabel('Cumulative fraction difference')
set(gca,'xtick',-1:0.5:1)
set(gca,'xticklabels',1:-0.5:-1)
axis square
legend('Same position','4 cm apart','8 cm apart','12 cm apart','16 cm apart','20 cm apart')
legend boxoff
set(gca,'fontsize',16)
box off
title('Figure 4H - CA1')

% Figure 4I - Difference in cumulative distributions for cells with different preferred positions (CA3):
figure
color_vec_CA3=[1 0.8 0; 1 0.4 0 ; 1 0 0; 0.8 0 0; 0.4 0 0 ; 0 0 0];
hold on
for n=1:6
    plot(corr_vec,conv(fliplr(mean_cdf_diff_diff_bins_within_day_pair_corr_CA3_novel(n,:)),smoothing_kernel,'same')./conv(ones(1,length(corr_vec)),smoothing_kernel,'same'),'-','linewidth',2,'color',color_vec_CA3(n,:))
end
plot([-1 0],[0 0],'--k','linewidth',2)
xlim([-1 0.5])
ylim([-0.05 0.1])
xlabel('Tuning curve correlation')
ylabel('Cumulative fraction difference')
set(gca,'xtick',-1:0.5:1)
set(gca,'xticklabels',1:-0.5:-1)
axis square
legend('Same position','4 cm apart','8 cm apart','12 cm apart','16 cm apart','20 cm apart')
legend boxoff
set(gca,'fontsize',16)
box off
title('Figure 4I - CA3')

%% Figure 5 - Place cells with peer-dependent tuning exhibit higher long-term tuning stability, but not higher tuning precision:

% Loading all the necessary data for Figure 5:
average_SI_SSR_bit_spike_dynamics_correlated_cells=nan(number_of_mice,maximal_number_of_sessions);
average_SI_SSR_bit_spike_dynamics_uncorrelated_cells=nan(number_of_mice,maximal_number_of_sessions);
PV_matrix_correlated_cells=nan(number_of_mice,maximal_number_of_sessions*num_trials,maximal_number_of_sessions*num_trials);
PV_matrix_uncorrelated_cells_control=nan(number_of_mice,maximal_number_of_sessions*num_trials,maximal_number_of_sessions*num_trials);
PV_dynamics_correlated_cells=nan(number_of_mice,maximal_number_of_sessions);
PV_dynamics_uncorrelated_cells_control=nan(number_of_mice,maximal_number_of_sessions);

within_day_tuning_corr_single_cells=cell(1,number_of_mice);
SI_SSR_single_cells=cell(1,number_of_mice);
correlation_percentile_single_cells=cell(1,number_of_mice);
across_days_SI_SSR_single_cells=cell(1,number_of_mice);
across_days_tuning_corr_single_cells=cell(1,number_of_mice);
across_days_correlation_percentile_single_cells=cell(1,number_of_mice);
for n=1:number_of_mice
    average_SI_SSR_bit_spike_dynamics_correlated_cells(n,:)=across_mice_data{n}.average_SI_SSR_bit_spike_dynamics_correlated_cells;
    average_SI_SSR_bit_spike_dynamics_uncorrelated_cells(n,:)=across_mice_data{n}.average_SI_SSR_bit_spike_dynamics_uncorrelated_cells;
    PV_matrix_correlated_cells(n,:,:)=across_mice_data{n}.PV_correlations_correlated_cells;
    PV_matrix_uncorrelated_cells_control(n,:,:)=across_mice_data{n}.PV_correlations_uncorrelated_cells_control;
    PV_dynamics_correlated_cells(n,:)=across_mice_data{n}.PV_dynamics_correlated_cells;
    PV_dynamics_uncorrelated_cells_control(n,:)=across_mice_data{n}.PV_dynamics_uncorrelated_cells_control;
    within_day_tuning_corr_single_cells{n}=across_mice_data{n}.within_day_tuning_corr_single_cells;
    SI_SSR_single_cells{n}=across_mice_data{n}.SI_SSR_single_cells;
    correlation_percentile_single_cells{n}=across_mice_data{n}.correlation_percentile_single_cells;
    across_days_tuning_corr_single_cells{n}=across_mice_data{n}.across_days_tuning_corr_single_cells;
    across_days_correlation_percentile_single_cells{n}=across_mice_data{n}.across_days_correlation_percentile_single_cells;
    across_days_SI_SSR_single_cells{n}=across_mice_data{n}.across_days_SI_SSR_single_cells;
end

% Figure 5A-C:
% Within day:
for n=1:number_of_mice
    within_day_tuning_corr_single_cells{n}(isnan(correlation_percentile_single_cells{n}))=[];
    SI_SSR_single_cells{n}(isnan(correlation_percentile_single_cells{n}))=[];
    correlation_percentile_single_cells{n}(isnan(correlation_percentile_single_cells{n}))=[];
    
    within_day_tuning_corr_single_cells{n}(isnan(SI_SSR_single_cells{n}))=[];
    correlation_percentile_single_cells{n}(isnan(SI_SSR_single_cells{n}))=[];
    SI_SSR_single_cells{n}(isnan(SI_SSR_single_cells{n}))=[];
    
    correlation_percentile_single_cells{n}(isnan(within_day_tuning_corr_single_cells{n}))=[];
    SI_SSR_single_cells{n}(isnan(within_day_tuning_corr_single_cells{n}))=[];
    within_day_tuning_corr_single_cells{n}(isnan(within_day_tuning_corr_single_cells{n}))=[];
end

information_vec=0.4:0.2:3;
average_info_vec=nan(1,length(information_vec)-1);
average_correlation_per_info=nan(number_of_mice,length(information_vec)-1);
for k=1:number_of_mice
    for n=1:length(information_vec)-1
        this_info_indexes=find(across_days_SI_SSR_single_cells{k}>=information_vec(n) & across_days_SI_SSR_single_cells{k}<=information_vec(n+1));
        average_correlation_per_info(k,n)=mean(across_days_tuning_corr_single_cells{k}(this_info_indexes),'omitnan');
        if k==1
            average_info_vec(n)=mean(information_vec(n:n+1));
        end
    end
end

mean_correlation_per_info_CA1=mean(average_correlation_per_info(mouse_group==1,:),'omitnan');
sd_correlation_per_info_CA1=std(average_correlation_per_info(mouse_group==1,:),'omitnan');
mean_correlation_per_info_CA3=mean(average_correlation_per_info(mouse_group==2,:),'omitnan');
sd_correlation_per_info_CA3=std(average_correlation_per_info(mouse_group==2,:),'omitnan');

dependence_vec=0:0.1:1;
average_dependence_vec=nan(1,length(dependence_vec)-1);
average_info_per_dependence=nan(number_of_mice,length(dependence_vec)-1);
average_within_correlation_per_dependence=nan(number_of_mice,length(dependence_vec)-1);
for k=1:number_of_mice
    for n=1:length(dependence_vec)-1
        this_dependence_indexes=find(correlation_percentile_single_cells{k}>=dependence_vec(n) & correlation_percentile_single_cells{k}<=dependence_vec(n+1));
        average_info_per_dependence(k,n)=mean(SI_SSR_single_cells{k}(this_dependence_indexes),'omitnan');
        average_within_correlation_per_dependence(k,n)=mean(within_day_tuning_corr_single_cells{k}(this_dependence_indexes),'omitnan');
        if k==1
            average_dependence_vec(n)=mean(dependence_vec(n:n+1));
        end
    end
end

mean_info_per_dependence_CA1=mean(average_info_per_dependence(mouse_group==1,:),'omitnan');
sd_info_per_dependence_CA1=std(average_info_per_dependence(mouse_group==1,:),'omitnan');
mean_info_per_dependence_CA3=mean(average_info_per_dependence(mouse_group==2,:),'omitnan');
sd_info_per_dependence_CA3=std(average_info_per_dependence(mouse_group==2,:),'omitnan');

mean_within_correlation_per_dependence_CA1=mean(average_within_correlation_per_dependence(mouse_group==1,:),'omitnan');
sd_within_correlation_per_dependence_CA1=std(average_within_correlation_per_dependence(mouse_group==1,:),'omitnan');
mean_within_correlation_per_dependence_CA3=mean(average_within_correlation_per_dependence(mouse_group==2,:),'omitnan');
sd_within_correlation_per_dependence_CA3=std(average_within_correlation_per_dependence(mouse_group==2,:),'omitnan');

% Across days:
for n=1:number_of_mice
    across_days_tuning_corr_single_cells{n}(isnan(across_days_SI_SSR_single_cells{n}))=[];
    across_days_correlation_percentile_single_cells{n}(isnan(across_days_SI_SSR_single_cells{n}))=[];
    across_days_SI_SSR_single_cells{n}(isnan(across_days_SI_SSR_single_cells{n}))=[];
    
    across_days_SI_SSR_single_cells{n}(isnan(across_days_tuning_corr_single_cells{n}))=[];
    across_days_correlation_percentile_single_cells{n}(isnan(across_days_tuning_corr_single_cells{n}))=[];
    across_days_tuning_corr_single_cells{n}(isnan(across_days_tuning_corr_single_cells{n}))=[];
    
    across_days_tuning_corr_single_cells{n}(isnan(across_days_correlation_percentile_single_cells{n}))=[];
    across_days_SI_SSR_single_cells{n}(isnan(across_days_correlation_percentile_single_cells{n}))=[];
    across_days_correlation_percentile_single_cells{n}(isnan(across_days_correlation_percentile_single_cells{n}))=[];
end

average_correlation_per_dependence=nan(number_of_mice,length(dependence_vec)-1);
for k=1:number_of_mice
    for n=1:length(dependence_vec)-1
        this_dependence_indexes=find(across_days_correlation_percentile_single_cells{k}>=dependence_vec(n) & across_days_correlation_percentile_single_cells{k}<=dependence_vec(n+1));
        average_correlation_per_dependence(k,n)=mean(across_days_tuning_corr_single_cells{k}(this_dependence_indexes),'omitnan');
    end
end
mean_correlation_per_dependence_CA1=mean(average_correlation_per_dependence(mouse_group==1,:),'omitnan');
sd_correlation_per_dependence_CA1=std(average_correlation_per_dependence(mouse_group==1,:),'omitnan');
mean_correlation_per_dependence_CA3=mean(average_correlation_per_dependence(mouse_group==2,:),'omitnan');
sd_correlation_per_dependence_CA3=std(average_correlation_per_dependence(mouse_group==2,:),'omitnan');

% Figure 5D and G - Spatial inforamtion for peer-dependent and
% peer-independent cells:
mean_SI_SSR_bit_spike_dynamics_correlated_cells_CA1=mean(average_SI_SSR_bit_spike_dynamics_correlated_cells(mouse_group==1,:),'omitnan');
std_SI_SSR_bit_spike_dynamics_correlated_cells_CA1=std(average_SI_SSR_bit_spike_dynamics_correlated_cells(mouse_group==1,:),'omitnan');
mean_SI_SSR_bit_spike_dynamics_correlated_cells_CA3=mean(average_SI_SSR_bit_spike_dynamics_correlated_cells(mouse_group==2,:),'omitnan');
std_SI_SSR_bit_spike_dynamics_correlated_cells_CA3=std(average_SI_SSR_bit_spike_dynamics_correlated_cells(mouse_group==2,:),'omitnan');

mean_SI_SSR_bit_spike_dynamics_uncorrelated_cells_CA1=mean(average_SI_SSR_bit_spike_dynamics_uncorrelated_cells(mouse_group==1,:),'omitnan');
std_SI_SSR_bit_spike_dynamics_uncorrelated_cells_CA1=std(average_SI_SSR_bit_spike_dynamics_uncorrelated_cells(mouse_group==1,:),'omitnan');
mean_SI_SSR_bit_spike_dynamics_uncorrelated_cells_CA3=mean(average_SI_SSR_bit_spike_dynamics_uncorrelated_cells(mouse_group==2,:),'omitnan');
std_SI_SSR_bit_spike_dynamics_uncorrelated_cells_CA3=std(average_SI_SSR_bit_spike_dynamics_uncorrelated_cells(mouse_group==2,:),'omitnan');

SI_SSR_per_mouse_correlated_cells=mean(average_SI_SSR_bit_spike_dynamics_correlated_cells,2,'omitnan');
SI_SSR_per_mouse_uncorrelated_cells=mean(average_SI_SSR_bit_spike_dynamics_uncorrelated_cells,2,'omitnan');

% Figure 5E and H - Within-day PV correlation for peer-dependent and
% peer-independent cells:
time_gap=2;
PV_within_day_correlated=nan(number_of_mice,maximal_number_of_sessions);
PV_within_day_uncorrelated_control=nan(number_of_mice,maximal_number_of_sessions);
for m=1:number_of_mice
    for n=1:maximal_number_of_sessions
        PV_within_day_correlated(m,n)=PV_matrix_correlated_cells(m,(n-1)*num_trials+1,(n-1)*num_trials+2);
        PV_within_day_uncorrelated_control(m,n)=PV_matrix_uncorrelated_cells_control(m,(n-1)*num_trials+1,(n-1)*num_trials+2);
    end
end

mean_PV_within_day_correlated_CA1=mean(PV_within_day_correlated(mouse_group==1,:),'omitnan');
std_PV_within_day_correlated_CA1=std(PV_within_day_correlated(mouse_group==1,:),'omitnan');
mean_PV_within_day_uncorrelated_CA1_control=mean(PV_within_day_uncorrelated_control(mouse_group==1,:),'omitnan');
std_PV_within_day_uncorrelated_CA1_control=std(PV_within_day_uncorrelated_control(mouse_group==1,:),'omitnan');
mean_PV_within_day_correlated_CA3=mean(PV_within_day_correlated(mouse_group==2,:),'omitnan');
std_PV_within_day_correlated_CA3=std(PV_within_day_correlated(mouse_group==2,:),'omitnan');
mean_PV_within_day_uncorrelated_CA3_control=mean(PV_within_day_uncorrelated_control(mouse_group==2,:),'omitnan');
std_PV_within_day_uncorrelated_CA3_control=std(PV_within_day_uncorrelated_control(mouse_group==2,:),'omitnan');

mean_per_mouse_PV_within_day_correlated=mean(PV_within_day_correlated,2,'omitnan');
mean_per_mouse_PV_within_day_uncorrelated_control=mean(PV_within_day_uncorrelated_control,2,'omitnan');

% Figure 5F and I - Across-days PV correlation for peer-dependent and
% peer-independent cells:
mean_PV_dynamics_correlated_cells_CA1=mean(PV_dynamics_correlated_cells(mouse_group==1,:),'omitnan');
std_PV_dynamics_correlated_cells_CA1=std(PV_dynamics_correlated_cells(mouse_group==1,:),'omitnan');
mean_PV_dynamics_uncorrelated_cells_CA1_control=mean(PV_dynamics_uncorrelated_cells_control(mouse_group==1,:),'omitnan');
std_PV_dynamics_uncorrelated_cells_CA1_control=std(PV_dynamics_uncorrelated_cells_control(mouse_group==1,:),'omitnan');
mean_PV_dynamics_correlated_cells_CA3=mean(PV_dynamics_correlated_cells(mouse_group==2,:),'omitnan');
std_PV_dynamics_correlated_cells_CA3=std(PV_dynamics_correlated_cells(mouse_group==2,:),'omitnan');
mean_PV_dynamics_uncorrelated_cells_CA3_control=mean(PV_dynamics_uncorrelated_cells_control(mouse_group==2,:),'omitnan');
std_PV_dynamics_uncorrelated_cells_CA3_control=std(PV_dynamics_uncorrelated_cells_control(mouse_group==2,:),'omitnan');

PV_consecutive_days_correlated=nan(number_of_mice,maximal_number_of_sessions);
PV_consecutive_days_uncorrelated_control=nan(number_of_mice,maximal_number_of_sessions);
for m=1:number_of_mice
    for n=1:maximal_number_of_sessions/2-1
        PV_consecutive_days_correlated(m,n+1)=mean([PV_matrix_correlated_cells(m,(n-1)*num_trials+1,n*num_trials+1), PV_matrix_correlated_cells(m,(n-1)*num_trials+2,n*num_trials+2)]);
        PV_consecutive_days_uncorrelated_control(m,n+1)=mean([PV_matrix_uncorrelated_cells_control(m,(n-1)*num_trials+1,n*num_trials+1), PV_matrix_uncorrelated_cells_control(m,(n-1)*num_trials+2,n*num_trials+2)]);
        PV_consecutive_days_correlated(m,n+maximal_number_of_sessions/2+1)=mean([PV_matrix_correlated_cells(m,(n+maximal_number_of_sessions/2-1)*num_trials+1,(n+maximal_number_of_sessions/2)*num_trials+1), PV_matrix_correlated_cells(m,(n+maximal_number_of_sessions/2-1)*num_trials+2,(n+maximal_number_of_sessions/2)*num_trials+2)]);
        PV_consecutive_days_uncorrelated_control(m,n+maximal_number_of_sessions/2+1)=mean([PV_matrix_uncorrelated_cells_control(m,(n+maximal_number_of_sessions/2-1)*num_trials+1,(n+maximal_number_of_sessions/2)*num_trials+1), PV_matrix_uncorrelated_cells_control(m,(n+maximal_number_of_sessions/2-1)*num_trials+2,(n+maximal_number_of_sessions/2)*num_trials+2)]);
    end
end
mean_per_mouse_PV_consecutive_days_correlated=mean(PV_consecutive_days_correlated,2,'omitnan');
mean_per_mouse_PV_consecutive_days_uncorrelated_control=mean(PV_consecutive_days_uncorrelated_control,2,'omitnan');

% Figure 5A - Tuning curve correlation as a function of spatial information:
figure
errorbar(average_info_vec,mean_correlation_per_info_CA1,sd_correlation_per_info_CA1./sqrt(sum(mouse_group==1)),'-*','markersize',5,'color','b','linewidth',2)
hold on
errorbar(average_info_vec,mean_correlation_per_info_CA3,sd_correlation_per_info_CA3./sqrt(sum(mouse_group==2)),'-*','markersize',5,'color','r','linewidth',2)
xlabel('Spatial information (bit/spike)')
ylabel('Tuning curve correlation')
ylim([0 1])
xlim([0.4 3])
box off
axis square
set(gca,'fontsize',16)
title('Figure 5A - Across days')

% Figure 5B - Spatial information as a function of tuning peer dependence:
figure
errorbar(average_dependence_vec,mean_info_per_dependence_CA1,sd_info_per_dependence_CA1./sqrt(sum(mouse_group==1)),'-*','markersize',5,'color','b','linewidth',2)
hold on
errorbar(average_dependence_vec,mean_info_per_dependence_CA3,sd_info_per_dependence_CA3./sqrt(sum(mouse_group==2)),'-*','markersize',5,'color','r','linewidth',2)
xlabel('Tuning peer-dependence')
ylabel('Spatial information (bit/spike)')
ylim([0 2])
box off
axis square
set(gca,'fontsize',16)
title('Figure 5B')

% Figure 5B, inset - Within-day tuning curve correlation as a function of
% tuning peer dependence:
figure
errorbar(average_dependence_vec,mean_within_correlation_per_dependence_CA1,sd_within_correlation_per_dependence_CA1./sqrt(sum(mouse_group==1)),'-*','markersize',5,'color','b','linewidth',2)
hold on
errorbar(average_dependence_vec,mean_within_correlation_per_dependence_CA3,sd_within_correlation_per_dependence_CA3./sqrt(sum(mouse_group==2)),'-*','markersize',5,'color','r','linewidth',2)
xlabel('Tuning peer-dependence')
ylabel('Tuning curve correlation')
ylim([0 1])
box off
axis square
set(gca,'fontsize',16)
title('Figure 5B, inset - Within day')

% Figure 5C - Across days tuning curve correlation as a function of
% tuning peer dependence:
figure
errorbar(average_dependence_vec,mean_correlation_per_dependence_CA1,sd_correlation_per_dependence_CA1./sqrt(sum(mouse_group==1)),'-*','markersize',5,'color','b','linewidth',2)
hold on
errorbar(average_dependence_vec,mean_correlation_per_dependence_CA3,sd_correlation_per_dependence_CA3./sqrt(sum(mouse_group==2)),'-*','markersize',5,'color','r','linewidth',2)
xlabel('Tuning peer-dependence')
ylabel('Tuning curve correlation')
ylim([0 1])
box off
axis square
legend('CA1','CA3')
legend boxoff
set(gca,'fontsize',16)
title('Figure 5C - Across days')

% Figure 5D - Spatial information over time for peer-dependent and
% peer-independent cells:
figure
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_SI_SSR_bit_spike_dynamics_correlated_cells_CA1(1:maximal_number_of_sessions/2),std_SI_SSR_bit_spike_dynamics_correlated_cells_CA1(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
hold on
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_SI_SSR_bit_spike_dynamics_uncorrelated_cells_CA1(1:maximal_number_of_sessions/2),std_SI_SSR_bit_spike_dynamics_uncorrelated_cells_CA1(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==1)),'--*b','linewidth',2)
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_SI_SSR_bit_spike_dynamics_correlated_cells_CA3(1:maximal_number_of_sessions/2),std_SI_SSR_bit_spike_dynamics_correlated_cells_CA3(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_SI_SSR_bit_spike_dynamics_uncorrelated_cells_CA3(1:maximal_number_of_sessions/2),std_SI_SSR_bit_spike_dynamics_uncorrelated_cells_CA3(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==2)),'--*r','linewidth',2)
errorbar(days_vec(maximal_number_of_sessions/2+1:end),mean_SI_SSR_bit_spike_dynamics_correlated_cells_CA1(maximal_number_of_sessions/2+1:end),std_SI_SSR_bit_spike_dynamics_correlated_cells_CA1(maximal_number_of_sessions/2+1:end)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
errorbar(days_vec(maximal_number_of_sessions/2+1:end),mean_SI_SSR_bit_spike_dynamics_uncorrelated_cells_CA1(maximal_number_of_sessions/2+1:end),std_SI_SSR_bit_spike_dynamics_uncorrelated_cells_CA1(maximal_number_of_sessions/2+1:end)./sqrt(sum(mouse_group==1)),'--*b','linewidth',2)
errorbar(days_vec(maximal_number_of_sessions/2+1:end),mean_SI_SSR_bit_spike_dynamics_correlated_cells_CA3(maximal_number_of_sessions/2+1:end),std_SI_SSR_bit_spike_dynamics_correlated_cells_CA3(maximal_number_of_sessions/2+1:end)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
errorbar(days_vec(maximal_number_of_sessions/2+1:end),mean_SI_SSR_bit_spike_dynamics_uncorrelated_cells_CA3(maximal_number_of_sessions/2+1:end),std_SI_SSR_bit_spike_dynamics_uncorrelated_cells_CA3(maximal_number_of_sessions/2+1:end)./sqrt(sum(mouse_group==2)),'--*r','linewidth',2)
plot([17 17],[0 800],'--','color','k','linewidth',2)
set(gca,'xtick',1:10:31)
xlim([1 33])
ylim([0 2])
xlabel('Time (days)')
ylabel('Spatial information (bit/spike)')
set(gca,'fontsize',16)
legend('Peer dependent - CA1','Peer independent - CA1','Peer dependent - CA3','Peer independent - CA3','Location','Southwest')
legend('boxoff')
box off
axis square
title('Figure 5D')

% Figure 5E - Within-day PV correlation over time for peer-dependent and
% peer-independent cells:
figure
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_PV_within_day_correlated_CA1(1:maximal_number_of_sessions/2),std_PV_within_day_correlated_CA1(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
hold on
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_PV_within_day_uncorrelated_CA1_control(1:maximal_number_of_sessions/2),std_PV_within_day_uncorrelated_CA1_control(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==1)),'--*b','linewidth',2)
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_PV_within_day_correlated_CA3(1:maximal_number_of_sessions/2),std_PV_within_day_correlated_CA3(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_PV_within_day_uncorrelated_CA3_control(1:maximal_number_of_sessions/2),std_PV_within_day_uncorrelated_CA3_control(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==2)),'--*r','linewidth',2)
errorbar(days_vec(maximal_number_of_sessions/2+1:end),mean_PV_within_day_correlated_CA1(maximal_number_of_sessions/2+1:end),std_PV_within_day_correlated_CA1(maximal_number_of_sessions/2+1:end)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
errorbar(days_vec(maximal_number_of_sessions/2+1:end),mean_PV_within_day_uncorrelated_CA1_control(maximal_number_of_sessions/2+1:end),std_PV_within_day_uncorrelated_CA1_control(maximal_number_of_sessions/2+1:end)./sqrt(sum(mouse_group==1)),'--*b','linewidth',2)
errorbar(days_vec(maximal_number_of_sessions/2+1:end),mean_PV_within_day_correlated_CA3(maximal_number_of_sessions/2+1:end),std_PV_within_day_correlated_CA3(maximal_number_of_sessions/2+1:end)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
errorbar(days_vec(maximal_number_of_sessions/2+1:end),mean_PV_within_day_uncorrelated_CA3_control(maximal_number_of_sessions/2+1:end),std_PV_within_day_uncorrelated_CA3_control(maximal_number_of_sessions/2+1:end)./sqrt(sum(mouse_group==2)),'--*r','linewidth',2)
plot([17 17],[0 800],'--','color','k','linewidth',2)
set(gca,'xtick',1:10:31)
xlim([1 33])
ylim([0 1])
xlabel('Time (days)')
ylabel('PV correlation')
set(gca,'fontsize',16)
legend('Peer dependent - CA1','Peer independent - CA1','Peer dependent - CA3','Peer independent - CA3','Location','Southwest')
legend('boxoff')
box off
axis square
title('Figure 5E - Within day')

% Figure 5F - PV correlation across days for peer-dependent and
% peer-independent cells:
figure
errorbar(2:2:14,mean_PV_dynamics_correlated_cells_CA1(2:8),std_PV_dynamics_correlated_cells_CA1(2:8)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
hold on
errorbar(2:2:14,mean_PV_dynamics_uncorrelated_cells_CA1_control(2:8),std_PV_dynamics_uncorrelated_cells_CA1_control(2:8)./sqrt(sum(mouse_group==1)),'--*b','linewidth',2)
errorbar(2:2:14,mean_PV_dynamics_correlated_cells_CA3(2:8),std_PV_dynamics_correlated_cells_CA3(2:8)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
errorbar(2:2:14,mean_PV_dynamics_uncorrelated_cells_CA3_control(2:8),std_PV_dynamics_uncorrelated_cells_CA3_control(2:8)./sqrt(sum(mouse_group==2)),'--*r','linewidth',2)
errorbar(0,mean_PV_dynamics_correlated_cells_CA1(1),std_PV_dynamics_correlated_cells_CA1(1)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
errorbar(0,mean_PV_dynamics_uncorrelated_cells_CA1_control(1),std_PV_dynamics_uncorrelated_cells_CA1_control(1)./sqrt(sum(mouse_group==1)),'--*b','linewidth',2)
errorbar(0,mean_PV_dynamics_correlated_cells_CA3(1),std_PV_dynamics_correlated_cells_CA3(1)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
errorbar(0,mean_PV_dynamics_uncorrelated_cells_CA3_control(1),std_PV_dynamics_uncorrelated_cells_CA3_control(1)./sqrt(sum(mouse_group==2)),'--*r','linewidth',2)
xlim([-0.5 14])
ylim([0 1])
xlabel('Elapsed time (days)')
ylabel('PV correlation')
legend('Peer dependent - CA1','Peer independent - CA1','Peer dependent - CA3','Peer independent - CA3','Location','Northeast')
legend('boxoff')
set(gca,'fontsize',16)
box off
axis square
title('Figure 5F - Across days')

% Figure 5G - Average spatial information for peer-dependent and
% peer-independent cells:
figure
x_vec_1=1+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_2=2+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_3=4+0.3*(0.5-rand(1,sum(mouse_group==2)));
x_vec_4=5+0.3*(0.5-rand(1,sum(mouse_group==2)));
plot([1 2],[mean(SI_SSR_per_mouse_correlated_cells(mouse_group==1)), mean(SI_SSR_per_mouse_uncorrelated_cells(mouse_group==1))],'.-b','markersize',30,'linewidth',2)
hold on
plot([4 5],[mean(SI_SSR_per_mouse_correlated_cells(mouse_group==2)), mean(SI_SSR_per_mouse_uncorrelated_cells(mouse_group==2))],'.-r','markersize',30,'linewidth',2)
errorbar([1 2],[mean(SI_SSR_per_mouse_correlated_cells(mouse_group==1)),mean(SI_SSR_per_mouse_uncorrelated_cells(mouse_group==1))],[std(SI_SSR_per_mouse_correlated_cells(mouse_group==1))./sqrt(sum(mouse_group==1)),std(SI_SSR_per_mouse_uncorrelated_cells(mouse_group==1))./sqrt(sum(mouse_group==1))],'.','linewidth',2,'color','b')
errorbar([4 5],[mean(SI_SSR_per_mouse_correlated_cells(mouse_group==2)),mean(SI_SSR_per_mouse_uncorrelated_cells(mouse_group==2))],[std(SI_SSR_per_mouse_correlated_cells(mouse_group==2))./sqrt(sum(mouse_group==2)),std(SI_SSR_per_mouse_uncorrelated_cells(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','r')
for n=1:sum(mouse_group==1)
    plot([x_vec_1(n), x_vec_2(n)],[SI_SSR_per_mouse_correlated_cells(CA1_indexes(n)),SI_SSR_per_mouse_uncorrelated_cells(CA1_indexes(n))],'-','color',[0.7 0.7 0.7])
end
for n=1:sum(mouse_group==2)
    plot([x_vec_3(n), x_vec_4(n)],[SI_SSR_per_mouse_correlated_cells(CA3_indexes(n)),SI_SSR_per_mouse_uncorrelated_cells(CA3_indexes(n))],'-','color',[0.7 0.7 0.7])
end
plot([1 2],[mean(SI_SSR_per_mouse_correlated_cells(mouse_group==1)), mean(SI_SSR_per_mouse_uncorrelated_cells(mouse_group==1))],'.-b','markersize',30,'linewidth',2)
plot([4 5],[mean(SI_SSR_per_mouse_correlated_cells(mouse_group==2)), mean(SI_SSR_per_mouse_uncorrelated_cells(mouse_group==2))],'.-r','markersize',30,'linewidth',2)
errorbar([1 2],[mean(SI_SSR_per_mouse_correlated_cells(mouse_group==1)),mean(SI_SSR_per_mouse_uncorrelated_cells(mouse_group==1))],[std(SI_SSR_per_mouse_correlated_cells(mouse_group==1))./sqrt(sum(mouse_group==1)),std(SI_SSR_per_mouse_uncorrelated_cells(mouse_group==1))./sqrt(sum(mouse_group==1))],'.','linewidth',2,'color','b')
errorbar([4 5],[mean(SI_SSR_per_mouse_correlated_cells(mouse_group==2)),mean(SI_SSR_per_mouse_uncorrelated_cells(mouse_group==2))],[std(SI_SSR_per_mouse_correlated_cells(mouse_group==2))./sqrt(sum(mouse_group==2)),std(SI_SSR_per_mouse_uncorrelated_cells(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','r')
xlim([0.5 5.5])
ylim([0 2])
ylabel('Spatial information (bit/spike)')
box off
hold on
axis square
set(gca,'xtick',[1 2 4 5])
set(gca,'xticklabels',{'Dep','Ind','Dep','Ind'})
set(gca,'fontsize',16)
axis square
box off
title('Figure 5G')

% Figure 5H - Average within-day PV correlation for peer-dependent and
% peer-independent cells:
figure
x_vec_1=1+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_2=2+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_3=4+0.3*(0.5-rand(1,sum(mouse_group==2)));
x_vec_4=5+0.3*(0.5-rand(1,sum(mouse_group==2)));
plot([1 2],[mean(mean_per_mouse_PV_within_day_correlated(mouse_group==1)), mean(mean_per_mouse_PV_within_day_uncorrelated_control(mouse_group==1))],'.-b','markersize',30,'linewidth',2)
hold on
plot([4 5],[mean(mean_per_mouse_PV_within_day_correlated(mouse_group==2)), mean(mean_per_mouse_PV_within_day_uncorrelated_control(mouse_group==2))],'.-r','markersize',30,'linewidth',2)
errorbar([1 2],[mean(mean_per_mouse_PV_within_day_correlated(mouse_group==1)),mean(mean_per_mouse_PV_within_day_uncorrelated_control(mouse_group==1))],[std(mean_per_mouse_PV_within_day_correlated(mouse_group==1))./sqrt(sum(mouse_group==1)),std(mean_per_mouse_PV_within_day_uncorrelated_control(mouse_group==1))./sqrt(sum(mouse_group==1))],'.','linewidth',2,'color','b')
errorbar([4 5],[mean(mean_per_mouse_PV_within_day_correlated(mouse_group==2)),mean(mean_per_mouse_PV_within_day_uncorrelated_control(mouse_group==2))],[std(mean_per_mouse_PV_within_day_correlated(mouse_group==2))./sqrt(sum(mouse_group==2)),std(mean_per_mouse_PV_within_day_uncorrelated_control(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','r')
for n=1:sum(mouse_group==1)
    plot([x_vec_1(n), x_vec_2(n)],[mean_per_mouse_PV_within_day_correlated(CA1_indexes(n)),mean_per_mouse_PV_within_day_uncorrelated_control(CA1_indexes(n))],'-','color',[0.7 0.7 0.7])
end
for n=1:sum(mouse_group==2)
    plot([x_vec_3(n), x_vec_4(n)],[mean_per_mouse_PV_within_day_correlated(CA3_indexes(n)),mean_per_mouse_PV_within_day_uncorrelated_control(CA3_indexes(n))],'-','color',[0.7 0.7 0.7])
end
plot([1 2],[mean(mean_per_mouse_PV_within_day_correlated(mouse_group==1)), mean(mean_per_mouse_PV_within_day_uncorrelated_control(mouse_group==1))],'.-b','markersize',30,'linewidth',2)
plot([4 5],[mean(mean_per_mouse_PV_within_day_correlated(mouse_group==2)), mean(mean_per_mouse_PV_within_day_uncorrelated_control(mouse_group==2))],'.-r','markersize',30,'linewidth',2)
errorbar([1 2],[mean(mean_per_mouse_PV_within_day_correlated(mouse_group==1)),mean(mean_per_mouse_PV_within_day_uncorrelated_control(mouse_group==1))],[std(mean_per_mouse_PV_within_day_correlated(mouse_group==1))./sqrt(sum(mouse_group==1)),std(mean_per_mouse_PV_within_day_uncorrelated_control(mouse_group==1))./sqrt(sum(mouse_group==1))],'.','linewidth',2,'color','b')
errorbar([4 5],[mean(mean_per_mouse_PV_within_day_correlated(mouse_group==2)),mean(mean_per_mouse_PV_within_day_uncorrelated_control(mouse_group==2))],[std(mean_per_mouse_PV_within_day_correlated(mouse_group==2))./sqrt(sum(mouse_group==2)),std(mean_per_mouse_PV_within_day_uncorrelated_control(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','r')
xlim([0.5 5.5])
ylim([0 1])
ylabel('PV correlation')
box off
hold on
axis square
set(gca,'xtick',[1 2 4 5])
set(gca,'xticklabels',{'Dep','Ind','Dep','Ind'})
set(gca,'fontsize',16)
axis square
box off
title('Figure 5H - Within day')

% Figure 5I - Average PV correlation across days for peer-dependent and
% peer-independent cells:
figure
x_vec_1=1+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_2=2+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_3=4+0.3*(0.5-rand(1,sum(mouse_group==2)));
x_vec_4=5+0.3*(0.5-rand(1,sum(mouse_group==2)));
plot([1 2],[mean(mean_per_mouse_PV_consecutive_days_correlated(mouse_group==1)), mean(mean_per_mouse_PV_consecutive_days_uncorrelated_control(mouse_group==1))],'.-b','markersize',30,'linewidth',2)
hold on
plot([4 5],[mean(mean_per_mouse_PV_consecutive_days_correlated(mouse_group==2)), mean(mean_per_mouse_PV_consecutive_days_uncorrelated_control(mouse_group==2))],'.-r','markersize',30,'linewidth',2)
errorbar([1 2],[mean(mean_per_mouse_PV_consecutive_days_correlated(mouse_group==1)),mean(mean_per_mouse_PV_consecutive_days_uncorrelated_control(mouse_group==1))],[std(mean_per_mouse_PV_consecutive_days_correlated(mouse_group==1))./sqrt(sum(mouse_group==1)),std(mean_per_mouse_PV_consecutive_days_uncorrelated_control(mouse_group==1))./sqrt(sum(mouse_group==1))],'.','linewidth',2,'color','b')
errorbar([4 5],[mean(mean_per_mouse_PV_consecutive_days_correlated(mouse_group==2)),mean(mean_per_mouse_PV_consecutive_days_uncorrelated_control(mouse_group==2))],[std(mean_per_mouse_PV_consecutive_days_correlated(mouse_group==2))./sqrt(sum(mouse_group==2)),std(mean_per_mouse_PV_consecutive_days_uncorrelated_control(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','r')
for n=1:sum(mouse_group==1)
    plot([x_vec_1(n), x_vec_2(n)],[mean_per_mouse_PV_consecutive_days_correlated(CA1_indexes(n)),mean_per_mouse_PV_consecutive_days_uncorrelated_control(CA1_indexes(n))],'-','color',[0.7 0.7 0.7])
end
for n=1:sum(mouse_group==2)
    plot([x_vec_3(n), x_vec_4(n)],[mean_per_mouse_PV_consecutive_days_correlated(CA3_indexes(n)),mean_per_mouse_PV_consecutive_days_uncorrelated_control(CA3_indexes(n))],'-','color',[0.7 0.7 0.7])
end
plot([1 2],[mean(mean_per_mouse_PV_consecutive_days_correlated(mouse_group==1)), mean(mean_per_mouse_PV_consecutive_days_uncorrelated_control(mouse_group==1))],'.-b','markersize',30,'linewidth',2)
plot([4 5],[mean(mean_per_mouse_PV_consecutive_days_correlated(mouse_group==2)), mean(mean_per_mouse_PV_consecutive_days_uncorrelated_control(mouse_group==2))],'.-r','markersize',30,'linewidth',2)
errorbar([1 2],[mean(mean_per_mouse_PV_consecutive_days_correlated(mouse_group==1)),mean(mean_per_mouse_PV_consecutive_days_uncorrelated_control(mouse_group==1))],[std(mean_per_mouse_PV_consecutive_days_correlated(mouse_group==1))./sqrt(sum(mouse_group==1)),std(mean_per_mouse_PV_consecutive_days_uncorrelated_control(mouse_group==1))./sqrt(sum(mouse_group==1))],'.','linewidth',2,'color','b')
errorbar([4 5],[mean(mean_per_mouse_PV_consecutive_days_correlated(mouse_group==2)),mean(mean_per_mouse_PV_consecutive_days_uncorrelated_control(mouse_group==2))],[std(mean_per_mouse_PV_consecutive_days_correlated(mouse_group==2))./sqrt(sum(mouse_group==2)),std(mean_per_mouse_PV_consecutive_days_uncorrelated_control(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','r')
xlim([0.5 5.5])
ylim([0 1])
ylabel('PV correlation')
box off
hold on
axis square
set(gca,'xtick',[1 2 4 5])
set(gca,'xticklabels',{'Dep','Ind','Dep','Ind'})
set(gca,'fontsize',16)
axis square
box off
title('Figure 5I - Across days')

%% Figure S1 - Cell detection and number of track traversals per day:

% Figure S1G-I:
number_of_neurons_dynamics=nan(number_of_mice,maximal_number_of_sessions);
average_firing_rates_dynamics=nan(number_of_mice,maximal_number_of_sessions);
number_of_track_traversals_dynamics=nan(number_of_mice,maximal_number_of_sessions);

for n=1:number_of_mice
    number_of_neurons_dynamics(n,:)=across_mice_data{n}.number_of_neurons_sessions;
    average_firing_rates_dynamics(n,:)=across_mice_data{n}.average_firing_rates_sessions;
    number_of_track_traversals_dynamics(n,:)=across_mice_data{n}.number_of_track_traversals_sessions;
end

average_number_of_neurons=mean(number_of_neurons_dynamics(:,:)');
average_firing_rates=mean(average_firing_rates_dynamics(:,:)');
average_number_of_track_traversals=mean(number_of_track_traversals_dynamics(:,:)');

% Figure S1G - Number of detected cells per day:
figure
x_vec_1=1+0.3*(0.5-rand(1,size(average_number_of_neurons(mouse_group==1),2)));
x_vec_2=2+0.3*(0.5-rand(1,size(average_number_of_neurons(mouse_group==2),2)));
bar([mean(average_number_of_neurons(mouse_group==1)),mean(average_number_of_neurons(mouse_group==2))],0.5,'FaceColor','none')
hold on
errorbar([mean(average_number_of_neurons(mouse_group==1)),mean(average_number_of_neurons(mouse_group==2))],[std(average_number_of_neurons(mouse_group==1))./sqrt(sum(mouse_group==1)),std(average_number_of_neurons(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',3,'color','k')
hold on
scatter(x_vec_1,average_number_of_neurons(mouse_group==1),25,'k','filled');
hold on
scatter(x_vec_2,average_number_of_neurons(mouse_group==2),25,'k','filled');
xlim([0.5 2.5])
ylim([0 800])
set(gca,'ytick',0:200:800)
ylabel('Number of detected cells')
set(gca,'fontsize',16)
box off
hold on
axis square
set(gca,'xtick',1:3)
set(gca,'xticklabels',{'CA1','CA3'})
set(gca,'fontsize',16)
set(gcf,'PaperPositionMode','auto')
axis square
box off
title('Figure S1G')

% Figure S1H - Average firing rates:
figure
x_vec_1=1+0.3*(0.5-rand(1,size(average_firing_rates(mouse_group==1),2)));
x_vec_2=2+0.3*(0.5-rand(1,size(average_firing_rates(mouse_group==2),2)));
bar([mean(average_firing_rates(mouse_group==1)),mean(average_firing_rates(mouse_group==2))],0.5,'FaceColor','none')
hold on
errorbar([mean(average_firing_rates(mouse_group==1)),mean(average_firing_rates(mouse_group==2))],[std(average_firing_rates(mouse_group==1))./sqrt(sum(mouse_group==1)),std(average_firing_rates(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',3,'color','k')
hold on
scatter(x_vec_1,average_firing_rates(mouse_group==1),25,'k','filled');
hold on
scatter(x_vec_2,average_firing_rates(mouse_group==2),25,'k','filled');
xlim([0.5 2.5])
ylim([0 0.5])
set(gca,'ytick',0:0.1:0.5)
ylabel('Average firing rate (spike/sec)')
box off
hold on
axis square
set(gca,'xtick',1:3)
set(gca,'xticklabels',{'CA1','CA3'})
set(gca,'fontsize',16)
set(gcf,'PaperPositionMode','auto')
axis square
box off
title('Figure S1H')

% Figure S1I - Number of track traversals:
figure
x_vec_1=1+0.3*(0.5-rand(1,size(average_number_of_track_traversals(mouse_group==1),2)));
x_vec_2=2+0.3*(0.5-rand(1,size(average_number_of_track_traversals(mouse_group==2),2)));
bar([mean(average_number_of_track_traversals(mouse_group==1)),mean(average_number_of_track_traversals(mouse_group==2))],0.5,'FaceColor','none')
hold on
errorbar([mean(average_number_of_track_traversals(mouse_group==1)),mean(average_number_of_track_traversals(mouse_group==2))],[std(average_number_of_track_traversals(mouse_group==1))./sqrt(sum(mouse_group==1)),std(average_number_of_track_traversals(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',3,'color','k')
hold on
scatter(x_vec_1,average_number_of_track_traversals(mouse_group==1),25,'k','filled');
hold on
scatter(x_vec_2,average_number_of_track_traversals(mouse_group==2),25,'k','filled');
xlim([0.5 2.5])
ylim([0 100])
set(gca,'ytick',0:20:100)
ylabel('Number of traversals')
box off
hold on
axis square
set(gca,'xtick',1:3)
set(gca,'xticklabels',{'CA1','CA3'})
set(gca,'fontsize',16)
set(gcf,'PaperPositionMode','auto')
axis square
box off
title('Figure S1I')

%% Figure S2 - Analyzing place field properties and recapitulating spatial code refinement
% by subsampling the data to obtain a fixed number of track traversals per day:

% Figure S2A-B - Place field properties:
within_field_threshold=0.3;
number_of_fields_dynamics=cell(number_of_mice,maximal_number_of_sessions,num_trials,2);
max_field_size_dynamics=cell(number_of_mice,maximal_number_of_sessions,num_trials,2);

for mouse_index=1:number_of_mice
    for n=1:maximal_number_of_sessions
        for k=1:num_trials
            this_trial_significant_cells_right=across_mice_data{mouse_index}.significant_right_trials{n,k};
            this_trial_rate_maps_right=squeeze(across_mice_data{mouse_index}.rate_maps_trials(n,k,this_trial_significant_cells_right,:,1));
            number_of_fields_dynamics{mouse_index,n,k,1}=nan(1,length(this_trial_significant_cells_right));
            max_field_size_dynamics{mouse_index,n,k,1}=nan(1,length(this_trial_significant_cells_right));
            if ~isempty(this_trial_significant_cells_right)
                for c=1:length(this_trial_significant_cells_right)
                    this_cell_rate_map=this_trial_rate_maps_right(c,:);
                    [max_rate_map,max_index_rate_map]=max(this_cell_rate_map);
                    this_cell_binarized_rate_map=this_cell_rate_map>within_field_threshold*max_rate_map;
                    this_cell_diff_rate_map=diff(this_cell_binarized_rate_map);
                    this_cell_num_fields=max([sum(this_cell_diff_rate_map==1) sum(this_cell_diff_rate_map==-1)]);
                    this_cell_total_field_size=sum(this_cell_binarized_rate_map);
                    if max_index_rate_map<20 && max_index_rate_map>1
                        end_of_main_field=find(this_cell_diff_rate_map(max_index_rate_map:end)==-1,1);
                        start_of_main_field=find(fliplr(this_cell_diff_rate_map(1:max_index_rate_map))==1,1);
                    elseif max_index_rate_map>1
                        end_of_main_field=2;
                        start_of_main_field=find(fliplr(this_cell_diff_rate_map(1:max_index_rate_map-1))==1,1);
                    elseif max_index_rate_map<20
                        end_of_main_field=find(this_cell_diff_rate_map(max_index_rate_map:end)==-1,1);
                        start_of_main_field=1;
                    end
                    if ~isempty(end_of_main_field) & ~isempty(start_of_main_field)
                        this_cell_main_field_size=end_of_main_field+start_of_main_field-2;
                    elseif ~isempty(end_of_main_field)
                        this_cell_main_field_size=end_of_main_field+max_index_rate_map-1;
                    elseif ~isempty(start_of_main_field)
                        this_cell_main_field_size=start_of_main_field+num_bins-max_index_rate_map-1;
                    else
                        this_cell_main_field_size=num_bins;
                    end
                    number_of_fields_dynamics{mouse_index,n,k,1}(c)=this_cell_num_fields;
                    max_field_size_dynamics{mouse_index,n,k,1}(c)=this_cell_main_field_size;
                end
            end
            
            this_trial_significant_cells_left=across_mice_data{mouse_index}.significant_left_trials{n,k};
            this_trial_rate_maps_left=squeeze(across_mice_data{mouse_index}.rate_maps_trials(n,k,this_trial_significant_cells_left,:,2));
            number_of_fields_dynamics{mouse_index,n,k,2}=nan(1,length(this_trial_significant_cells_left));
            max_field_size_dynamics{mouse_index,n,k,2}=nan(1,length(this_trial_significant_cells_left));
            if ~isempty(this_trial_significant_cells_left)
                for c=1:length(this_trial_significant_cells_left)
                    this_cell_rate_map=this_trial_rate_maps_left(c,:);
                    [max_rate_map,max_index_rate_map]=max(this_cell_rate_map);
                    this_cell_binarized_rate_map=this_cell_rate_map>within_field_threshold*max_rate_map;
                    this_cell_diff_rate_map=diff(this_cell_binarized_rate_map);
                    this_cell_num_fields=max([sum(this_cell_diff_rate_map==1) sum(this_cell_diff_rate_map==-1)]);
                    this_cell_total_field_size=sum(this_cell_binarized_rate_map);
                    if max_index_rate_map<20 && max_index_rate_map>1
                        end_of_main_field=find(this_cell_diff_rate_map(max_index_rate_map:end)==-1,1);
                        start_of_main_field=find(fliplr(this_cell_diff_rate_map(1:max_index_rate_map))==1,1);
                    elseif max_index_rate_map>1
                        end_of_main_field=2;
                        start_of_main_field=find(fliplr(this_cell_diff_rate_map(1:max_index_rate_map-1))==1,1);
                    elseif max_index_rate_map<20
                        end_of_main_field=find(this_cell_diff_rate_map(max_index_rate_map:end)==-1,1);
                        start_of_main_field=1;
                    end
                    if ~isempty(end_of_main_field) & ~isempty(start_of_main_field)
                        this_cell_main_field_size=end_of_main_field+start_of_main_field-2;
                    elseif ~isempty(end_of_main_field)
                        this_cell_main_field_size=end_of_main_field+max_index_rate_map-1;
                    elseif ~isempty(start_of_main_field)
                        this_cell_main_field_size=start_of_main_field+num_bins-max_index_rate_map-1;
                    else
                        this_cell_main_field_size=num_bins;
                    end
                    number_of_fields_dynamics{mouse_index,n,k,2}(c)=this_cell_num_fields;
                    max_field_size_dynamics{mouse_index,n,k,2}(c)=this_cell_main_field_size;
                end
            end
        end
    end
end

average_number_of_fields_dynamics=nan(number_of_mice,maximal_number_of_sessions,num_trials);
average_max_field_size_dynamics=nan(number_of_mice,maximal_number_of_sessions,num_trials);
for mouse_index=1:number_of_mice
    for n=1:maximal_number_of_sessions
        for k=1:num_trials
            average_number_of_fields_dynamics(mouse_index,n,k)=mean([number_of_fields_dynamics{mouse_index,n,k,1},number_of_fields_dynamics{mouse_index,n,k,2}],'omitnan');
            average_max_field_size_dynamics(mouse_index,n,k)=mean([max_field_size_dynamics{mouse_index,n,k,1},max_field_size_dynamics{mouse_index,n,k,2}],'omitnan');
        end
    end
end
average_number_of_fields_days=mean(average_number_of_fields_dynamics,3,'omitnan');
average_max_field_size_days=mean(average_max_field_size_dynamics,3,'omitnan');

average_number_of_fields_novel=mean(average_number_of_fields_days(:,novel_sessions),2,'omitnan');
average_max_field_size_novel=mean(average_max_field_size_days(:,novel_sessions),2,'omitnan');

% Figure S2C-E - Comparison of the spatial code across environments:
SI_SSR_bit_spike_dynamics_A=SI_SSR_bit_spike_dynamics(:,1:maximal_number_of_sessions/2);
SI_SSR_bit_spike_dynamics_B=SI_SSR_bit_spike_dynamics(:,maximal_number_of_sessions/2+1:end);
SI_SSR_bit_spike_dynamics_diff=SI_SSR_bit_spike_dynamics_B-SI_SSR_bit_spike_dynamics_A;
mean_SI_SSR_bit_spike_dynamics_diff_CA1=mean(SI_SSR_bit_spike_dynamics_diff(mouse_group==1,:),'omitnan');
std_SI_SSR_bit_spike_dynamics_diff_CA1=std(SI_SSR_bit_spike_dynamics_diff(mouse_group==1,:),'omitnan');
mean_SI_SSR_bit_spike_dynamics_diff_CA3=mean(SI_SSR_bit_spike_dynamics_diff(mouse_group==2,:),'omitnan');
std_SI_SSR_bit_spike_dynamics_diff_CA3=std(SI_SSR_bit_spike_dynamics_diff(mouse_group==2,:),'omitnan');

across_sides_PV_dynamics_A=PV_correlations_across_sides_maximal_shift(:,1:maximal_number_of_sessions/2);
across_sides_PV_dynamics_B=PV_correlations_across_sides_maximal_shift(:,maximal_number_of_sessions/2+1:end);
across_sides_PV_dynamics_diff=across_sides_PV_dynamics_B-across_sides_PV_dynamics_A;
mean_across_sides_PV_dynamics_diff_CA1=mean(across_sides_PV_dynamics_diff(mouse_group==1,:),'omitnan');
std_across_sides_PV_dynamics_diff_CA1=std(across_sides_PV_dynamics_diff(mouse_group==1,:),'omitnan');
mean_across_sides_PV_dynamics_diff_CA3=mean(across_sides_PV_dynamics_diff(mouse_group==2,:),'omitnan');
std_across_sides_PV_dynamics_diff_CA3=std(across_sides_PV_dynamics_diff(mouse_group==2,:),'omitnan');

within_day_PV_dynamics_A=within_day_PV_correlations(:,1:maximal_number_of_sessions/2);
within_day_PV_dynamics_B=within_day_PV_correlations(:,maximal_number_of_sessions/2+1:end);
within_day_PV_dynamics_diff=within_day_PV_dynamics_B-within_day_PV_dynamics_A;
mean_within_day_PV_dynamics_diff_CA1=mean(within_day_PV_dynamics_diff(mouse_group==1,:),'omitnan');
std_within_day_PV_dynamics_diff_CA1=std(within_day_PV_dynamics_diff(mouse_group==1,:),'omitnan');
mean_within_day_PV_dynamics_diff_CA3=mean(within_day_PV_dynamics_diff(mouse_group==2,:),'omitnan');
std_within_day_PV_dynamics_diff_CA3=std(within_day_PV_dynamics_diff(mouse_group==2,:),'omitnan');

% Figure S2A - Average place field size:
figure
x_vec_1=1+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_2=2+0.3*(0.5-rand(1,sum(mouse_group==2)));
bar(bin_size*[mean(average_max_field_size_novel(mouse_group==1)),mean(average_max_field_size_novel(mouse_group==2))],0.5,'FaceColor','none')
hold on
errorbar(bin_size*[mean(average_max_field_size_novel(mouse_group==1)),mean(average_max_field_size_novel(mouse_group==2))],bin_size*[std(average_max_field_size_novel(mouse_group==1))./sqrt(sum(mouse_group==1)),std(average_max_field_size_novel(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',3,'color','k')
hold on
scatter(x_vec_1,bin_size*average_max_field_size_novel(mouse_group==1),25,'k','filled');
hold on
scatter(x_vec_2,bin_size*average_max_field_size_novel(mouse_group==2),25,'k','filled');
xlim([0.5 2.5])
ylim([0 15])
set(gca,'ytick',0:5:20)
ylabel('Place field size (cm)')
set(gca,'fontsize',16)
box off
hold on
axis square
set(gca,'xticklabels',{'CA1','CA3'})
set(gca,'fontsize',16)
set(gcf,'PaperPositionMode','auto')
axis square
box off
title('Figure S2A')

% Figure S2B - Average number of place fields per cell:
figure
x_vec_1=1+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_2=2+0.3*(0.5-rand(1,sum(mouse_group==2)));
bar([mean(average_number_of_fields_novel(mouse_group==1)),mean(average_number_of_fields_novel(mouse_group==2))],0.5,'FaceColor','none')
hold on
errorbar([mean(average_number_of_fields_novel(mouse_group==1)),mean(average_number_of_fields_novel(mouse_group==2))],[std(average_number_of_fields_novel(mouse_group==1))./sqrt(sum(mouse_group==1)),std(average_number_of_fields_novel(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',3,'color','k')
hold on
scatter(x_vec_1,average_number_of_fields_novel(mouse_group==1),25,'k','filled');
hold on
scatter(x_vec_2,average_number_of_fields_novel(mouse_group==2),25,'k','filled');
xlim([0.5 2.5])
ylim([0 3])
set(gca,'ytick',0:3)
ylabel('Number of place fields per cell')
set(gca,'fontsize',16)
box off
hold on
axis square
set(gca,'xticklabels',{'CA1','CA3'})
set(gca,'fontsize',16)
set(gcf,'PaperPositionMode','auto')
axis square
box off
title('Figure S2B')

% Figure S2C - Comparison of spatial information across environments:
figure
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_SI_SSR_bit_spike_dynamics_diff_CA1(1:maximal_number_of_sessions/2),std_SI_SSR_bit_spike_dynamics_diff_CA1(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
hold on
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_SI_SSR_bit_spike_dynamics_diff_CA3(1:maximal_number_of_sessions/2),std_SI_SSR_bit_spike_dynamics_diff_CA3(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
plot([1 15], [0 0],'--k','linewidth',2)
xlim([1 15])
ylim([-0.5 0.5])
xlabel('Time (days)')
ylabel('\delta Spatial information (bit/spike)')
set(gca,'fontsize',16)
legend('CA1','CA3')
legend('boxoff')
box off
axis square
title('Figure S2C - B minus A')

% Figure S2D - Comparison of across directions PV correlation across environments:
figure
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_across_sides_PV_dynamics_diff_CA1(1:maximal_number_of_sessions/2),std_across_sides_PV_dynamics_diff_CA1(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
hold on
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_across_sides_PV_dynamics_diff_CA3(1:maximal_number_of_sessions/2),std_across_sides_PV_dynamics_diff_CA3(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
plot([1 15], [0 0],'--k','linewidth',2)
xlim([1 15])
ylim([-0.5 0.5])
xlabel('Time (days)')
ylabel('\delta PV correlation across directions')
set(gca,'fontsize',16)
legend('CA1','CA3')
legend('boxoff')
box off
axis square
title('Figure S2D - B minus A')

% Figure S2E - Comparison of within-day PV correlation across environments:
figure
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_within_day_PV_dynamics_diff_CA1(1:maximal_number_of_sessions/2),std_within_day_PV_dynamics_diff_CA1(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
hold on
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_within_day_PV_dynamics_diff_CA3(1:maximal_number_of_sessions/2),std_within_day_PV_dynamics_diff_CA3(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
plot([1 15], [0 0],'--k','linewidth',2)
xlim([1 15])
ylim([-0.5 0.5])
xlabel('Time (days)')
ylabel('\delta PV correlation within day')
set(gca,'fontsize',16)
legend('CA1','CA3')
legend('boxoff')
box off
axis square
title('Figure S2E - B minus A')

% Figure S2F-I - recapitulating spatial code refinement by subsampling
% the data to obtain a fixed number of track traversals per day:
% Figure S2F-G - Spatial information:
SI_SSR_bit_spike_place_cells_novel_CA1=[];
SI_SSR_bit_spike_place_cells_novel_CA3=[];
for n=1:number_of_mice
    for k=1:maximal_number_of_sessions
        if mouse_group(n)==1
            if ~isempty(intersect(k,novel_sessions))
                SI_SSR_bit_spike_place_cells_novel_CA1=[SI_SSR_bit_spike_place_cells_novel_CA1,across_mice_data{n}.SI_SSR_bit_spike_sessions_subsampled{k}];
            end
        elseif mouse_group(n)==2
            if ~isempty(intersect(k,novel_sessions))
                SI_SSR_bit_spike_place_cells_novel_CA3=[SI_SSR_bit_spike_place_cells_novel_CA3,across_mice_data{n}.SI_SSR_bit_spike_sessions_subsampled{k}];
            end
        end
    end
end

SI_SSR_bit_spike_dynamics=nan(number_of_mice,maximal_number_of_sessions);
for n=1:number_of_mice
    SI_SSR_bit_spike_dynamics(n,:)=across_mice_data{n}.average_SI_SSR_bit_spike_sessions_subsampled;
end

mean_SI_SSR_bit_spike_dynamics_CA1=mean(SI_SSR_bit_spike_dynamics(mouse_group==1,:),'omitnan');
std_SI_SSR_bit_spike_dynamics_CA1=std(SI_SSR_bit_spike_dynamics(mouse_group==1,:),'omitnan');
mean_SI_SSR_bit_spike_dynamics_CA3=mean(SI_SSR_bit_spike_dynamics(mouse_group==2,:),'omitnan');
std_SI_SSR_bit_spike_dynamics_CA3=std(SI_SSR_bit_spike_dynamics(mouse_group==2,:),'omitnan');
average_SI_SSR_bit_spike=mean(SI_SSR_bit_spike_dynamics(:,novel_sessions)');

% Figure S2F - Distribution of spatial information:
figure
x_vec=-0.5:0.2:4.9;
[n1,~]=hist(SI_SSR_bit_spike_place_cells_novel_CA1,x_vec);
[n2,~]=hist(SI_SSR_bit_spike_place_cells_novel_CA3,x_vec);
n1=n1./sum(n1);
n2=n2./sum(n2);
plot(x_vec,n1,'-b','linewidth',2);
hold on
plot(x_vec,n2,'-r','linewidth',2);
xlim([0 4])
xlabel('Spatial information (bit/spike)')
ylabel('Fraction of cells')
legend('CA1','CA3')
legend boxoff
set(gca,'fontsize',16)
box off
axis square
title('Figure S2F')

% Figure S2F, inset - Average spatial information:
figure
x_vec_1=1+0.3*(0.5-rand(1,size(average_SI_SSR_bit_spike(mouse_group==1),2)));
x_vec_2=2+0.3*(0.5-rand(1,size(average_SI_SSR_bit_spike(mouse_group==2),2)));
bar([mean(average_SI_SSR_bit_spike(mouse_group==1)),mean(average_SI_SSR_bit_spike(mouse_group==2))],0.5,'FaceColor','none')
hold on
errorbar([mean(average_SI_SSR_bit_spike(mouse_group==1)),mean(average_SI_SSR_bit_spike(mouse_group==2))],[std(average_SI_SSR_bit_spike(mouse_group==1))./sqrt(sum(mouse_group==1)),std(average_SI_SSR_bit_spike(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',3,'color','k')
hold on
scatter(x_vec_1,average_SI_SSR_bit_spike(mouse_group==1),15,'k','filled');
hold on
scatter(x_vec_2,average_SI_SSR_bit_spike(mouse_group==2),15,'k','filled');
xlim([0.5 2.5])
ylim([0 2])
set(gca,'ytick',0.5:0.5:2)
ylabel('Spatial information (bit/spike)')
set(gca,'fontsize',16)
box off
hold on
axis square
set(gca,'xtick',1:3)
set(gca,'xticklabels',{'CA1','CA3'})
title('Figure S2F, inset')
set(gca,'fontsize',16)
set(gcf,'PaperPositionMode','auto')
axis square
box off

% Figure S2G - Spatial information over days:
figure
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_SI_SSR_bit_spike_dynamics_CA1(1:maximal_number_of_sessions/2),std_SI_SSR_bit_spike_dynamics_CA1(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
hold on
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_SI_SSR_bit_spike_dynamics_CA3(1:maximal_number_of_sessions/2),std_SI_SSR_bit_spike_dynamics_CA3(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
hold on
errorbar(days_vec(maximal_number_of_sessions/2+1:end),mean_SI_SSR_bit_spike_dynamics_CA1(maximal_number_of_sessions/2+1:end),std_SI_SSR_bit_spike_dynamics_CA1(maximal_number_of_sessions/2+1:end)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
hold on
errorbar(days_vec(maximal_number_of_sessions/2+1:end),mean_SI_SSR_bit_spike_dynamics_CA3(maximal_number_of_sessions/2+1:end),std_SI_SSR_bit_spike_dynamics_CA3(maximal_number_of_sessions/2+1:end)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
plot([17 17],[0 800],'--','color','k','linewidth',2)
set(gca,'xtick',1:10:31)
xlim([1 33])
ylim([0 2])
xlabel('Time (days)')
ylabel('Spatial information (bit/spike)')
title('Figure S2G')
set(gca,'fontsize',16)
legend('CA1','CA3','Location','Southwest')
legend('boxoff')
box off
axis square

% Figure S2H-I - PV correlations:
per_day_PV_correlations_bin_resolution=nan(number_of_mice,maximal_number_of_sessions,2*num_bins*2,2*num_bins*2);
for k=1:number_of_mice
    this_mouse_PV_correlations_bin_resolution=across_mice_data{k}.PV_correlations_bin_resolution_subsampled;
    for n=1:maximal_number_of_sessions
        per_day_PV_correlations_bin_resolution(k,n,:,:)=this_mouse_PV_correlations_bin_resolution((n-1)*2*num_bins*2+1:n*2*num_bins*2,(n-1).*2*num_bins*2+1:n.*2*num_bins*2);
    end
end

per_day_PV_correlations_bin_resolution_CA1=squeeze(mean(per_day_PV_correlations_bin_resolution(mouse_group==1,:,:,:),'omitnan'));
per_day_PV_correlations_bin_resolution_CA3=squeeze(mean(per_day_PV_correlations_bin_resolution(mouse_group==2,:,:,:),'omitnan'));
PV_correlations_bin_resolution_novel_CA1=squeeze(mean(per_day_PV_correlations_bin_resolution_CA1(novel_sessions,:,:),'omitnan'));
PV_correlations_bin_resolution_novel_CA3=squeeze(mean(per_day_PV_correlations_bin_resolution_CA3(novel_sessions,:,:),'omitnan'));

within_day_PV_correlations=nan(number_of_mice,maximal_number_of_sessions);
for n=1:number_of_mice
    within_day_PV_correlations(n,:)=across_mice_data{n}.within_day_PV_correlations_subsampled;
end
mean_within_day_PV_correlations_CA1=mean(within_day_PV_correlations(mouse_group==1,:),'omitnan');
std_within_day_PV_correlations_CA1=std(within_day_PV_correlations(mouse_group==1,:),'omitnan');
mean_within_day_PV_correlations_CA3=mean(within_day_PV_correlations(mouse_group==2,:),'omitnan');
std_within_day_PV_correlations_CA3=std(within_day_PV_correlations(mouse_group==2,:),'omitnan');

average_within_day_PV_correlations=mean(within_day_PV_correlations(:,novel_sessions)');

% PV correlation across directions:
PV_correlations_across_sides_maximal_shift=nan(number_of_mice,maximal_number_of_sessions);
possible_shifts=-5:5;
for n=1:number_of_mice
    for k=1:maximal_number_of_sessions
        this_day_PV_correlations=squeeze(per_day_PV_correlations_bin_resolution(n,k,:,:));
        
        temp_PV_correlations_vs_shift=nan(1,length(possible_shifts));
        for shift_ind=1:length(possible_shifts)
            temp_PV_correlations_this_shift=nan(1,2*num_bins);
            for s=1:num_bins-abs(possible_shifts(shift_ind))
                if possible_shifts(shift_ind)>=0
                    temp_PV_correlations_this_shift(s)=this_day_PV_correlations(s,s+num_bins+possible_shifts(shift_ind));
                    temp_PV_correlations_this_shift(num_bins+s)=this_day_PV_correlations(2*num_bins+s,2*num_bins+s+num_bins+possible_shifts(shift_ind));
                else
                    temp_PV_correlations_this_shift(s)=this_day_PV_correlations(s-possible_shifts(shift_ind),s+num_bins);
                    temp_PV_correlations_this_shift(num_bins+s)=this_day_PV_correlations(2*num_bins+s-possible_shifts(shift_ind),2*num_bins+s+num_bins);
                end
            end
            temp_PV_correlations_vs_shift(shift_ind)=mean(temp_PV_correlations_this_shift,'omitnan');
        end
        [this_maximal_shift_PV_correlation,~]=max(temp_PV_correlations_vs_shift);
        PV_correlations_across_sides_maximal_shift(n,k)=this_maximal_shift_PV_correlation;
    end
end

PV_correlations_novel_maximal_shift=mean(PV_correlations_across_sides_maximal_shift(:,novel_sessions),2,'omitnan');

mean_across_sides_maximal_shift_PV_dynamics_CA1=mean(PV_correlations_across_sides_maximal_shift(mouse_group==1,:),'omitnan');
std_across_sides_maximal_shift_PV_dynamics_CA1=std(PV_correlations_across_sides_maximal_shift(mouse_group==1,:),'omitnan');
mean_across_sides_maximal_shift_PV_dynamics_CA3=mean(PV_correlations_across_sides_maximal_shift(mouse_group==2,:),'omitnan');
std_across_sides_maximal_shift_PV_dynamics_CA3=std(PV_correlations_across_sides_maximal_shift(mouse_group==2,:),'omitnan');

PV_correlations_novel_maximal_shift_L=nan(1,number_of_mice);
PV_correlations_novel_maximal_shift_straight=nan(1,number_of_mice);
for n=1:number_of_mice
    if n>=3 && n<=5
        PV_correlations_novel_maximal_shift_L(n)=mean(PV_correlations_across_sides_maximal_shift(n,1:3),2,'omitnan');
        PV_correlations_novel_maximal_shift_straight(n)=mean(PV_correlations_across_sides_maximal_shift(n,9:11),2,'omitnan');
    else
        PV_correlations_novel_maximal_shift_L(n)=mean(PV_correlations_across_sides_maximal_shift(n,9:11),2,'omitnan');
        PV_correlations_novel_maximal_shift_straight(n)=mean(PV_correlations_across_sides_maximal_shift(n,1:3),2,'omitnan');
    end
end

% Figure S2H - PV correlation across directions over time:
figure
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_across_sides_maximal_shift_PV_dynamics_CA1(1:maximal_number_of_sessions/2),std_across_sides_maximal_shift_PV_dynamics_CA1(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
hold on
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_across_sides_maximal_shift_PV_dynamics_CA3(1:maximal_number_of_sessions/2),std_across_sides_maximal_shift_PV_dynamics_CA3(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
errorbar(days_vec(maximal_number_of_sessions/2+1:end),mean_across_sides_maximal_shift_PV_dynamics_CA1(maximal_number_of_sessions/2+1:end),std_across_sides_maximal_shift_PV_dynamics_CA1(maximal_number_of_sessions/2+1:end)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
errorbar(days_vec(maximal_number_of_sessions/2+1:end),mean_across_sides_maximal_shift_PV_dynamics_CA3(maximal_number_of_sessions/2+1:end),std_across_sides_maximal_shift_PV_dynamics_CA3(maximal_number_of_sessions/2+1:end)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
plot([17 17],[-0.1 1],'--','color','k','linewidth',2)
plot([1 33],[0 0],'--','color','k','linewidth',2)
set(gca,'xtick',1:10:31)
xlim([1 33])
ylim([-0.1 0.3])
set(gca,'ytick',-0.1:0.1:0.3)
xlabel('Time (days)')
ylabel('PV ocrrelations across directions')
set(gca,'fontsize',16)
legend('CA1','CA3')
legend boxoff
box off
axis square
title('Figure S2H - Within session')

% Figure S2I - Within-day PV correlations over time:
figure
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_within_day_PV_correlations_CA1(1:maximal_number_of_sessions/2),std_within_day_PV_correlations_CA1(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
hold on
errorbar(days_vec(1:maximal_number_of_sessions/2),mean_within_day_PV_correlations_CA3(1:maximal_number_of_sessions/2),std_within_day_PV_correlations_CA3(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
errorbar(days_vec(maximal_number_of_sessions/2+1:end),mean_within_day_PV_correlations_CA1(maximal_number_of_sessions/2+1:end),std_within_day_PV_correlations_CA1(maximal_number_of_sessions/2+1:end)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
errorbar(days_vec(maximal_number_of_sessions/2+1:end),mean_within_day_PV_correlations_CA3(maximal_number_of_sessions/2+1:end),std_within_day_PV_correlations_CA3(maximal_number_of_sessions/2+1:end)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
plot([17 17],[0 800],'--','color','k','linewidth',2)
set(gca,'xtick',1:10:31)
xlim([1 33])
ylim([0 1])
set(gca,'ytick',-0.2:0.2:1)
xlabel('Time (days)')
ylabel('PV correlation')
set(gca,'fontsize',16)
box off
axis square
legend('CA1','CA3','Location','southeast')
legend('boxoff')
title('Figure S2I - Within day')

%% Figure S3 - Tracking the same neurons across multiple days of imaging:

% Loading the data necessary for Figure S3:
number_of_registered_neurons=nan(1,number_of_mice);
number_of_neurons_dynamics=nan(number_of_mice,maximal_number_of_sessions);
overall_registration_errors=nan(1,number_of_mice);
PV_dynamics_subsampled=nan(number_of_mice,maximal_number_of_sessions);
rate_dynamics_subsampled=nan(number_of_mice,maximal_number_of_sessions);
for n=1:number_of_mice
    number_of_neurons_dynamics(n,:)=across_mice_data{n}.number_of_neurons_sessions;
    number_of_registered_neurons(n)=size(across_mice_data{n}.cell_to_index_map,1);
    overall_registration_errors(n)=across_mice_data{n}.overall_registration_errors;
    PV_dynamics_subsampled(n,:)=across_mice_data{n}.PV_dynamics_subsampled(1:maximal_number_of_sessions);
    rate_dynamics_subsampled(n,:)=across_mice_data{n}.rate_dynamics_subsampled(1:maximal_number_of_sessions);
end
average_number_of_neurons=mean(number_of_neurons_dynamics(:,:)');
fraction_of_neurons=(average_number_of_neurons./number_of_registered_neurons)*100;
overall_registration_errors_CA1=overall_registration_errors(mouse_group==1);
overall_registration_errors_CA3=overall_registration_errors(mouse_group==2);

fraction_of_neurons_dynamics=number_of_neurons_dynamics./number_of_registered_neurons';
mean_fraction_of_neurons_dynamics_CA1=mean(fraction_of_neurons_dynamics(mouse_group==1,:),'omitnan');
std_fraction_of_neurons_dynamics_CA1=std(fraction_of_neurons_dynamics(mouse_group==1,:),'omitnan');
mean_fraction_of_neurons_dynamics_CA3=mean(fraction_of_neurons_dynamics(mouse_group==2,:),'omitnan');
std_fraction_of_neurons_dynamics_CA3=std(fraction_of_neurons_dynamics(mouse_group==2,:),'omitnan');

mean_PV_dynamics_CA1=mean(PV_dynamics_subsampled(mouse_group==1,:),'omitnan');
std_PV_dynamics_CA1=std(PV_dynamics_subsampled(mouse_group==1,:),'omitnan');
mean_rate_dynamics_CA1=mean(rate_dynamics_subsampled(mouse_group==1,:),'omitnan');
std_rate_dynamics_CA1=std(rate_dynamics_subsampled(mouse_group==1,:),'omitnan');

mean_PV_dynamics_CA3=mean(PV_dynamics(mouse_group==2,:),'omitnan');
std_PV_dynamics_CA3=std(PV_dynamics(mouse_group==2,:),'omitnan');
mean_rate_dynamics_CA3=mean(rate_dynamics_subsampled(mouse_group==2,:),'omitnan');
std_rate_dynamics_CA3=std(rate_dynamics_subsampled(mouse_group==2,:),'omitnan');

% Figure S3D-F - Cell registration quality:
% Figure S3D - Total number of registered neurons per mouse:
figure
x_vec_1=1+0.3*(0.5-rand(1,size(number_of_registered_neurons(mouse_group==1),2)));
x_vec_2=2+0.3*(0.5-rand(1,size(number_of_registered_neurons(mouse_group==2),2)));
bar([mean(number_of_registered_neurons(mouse_group==1)),mean(number_of_registered_neurons(mouse_group==2))],0.5,'FaceColor','none')
hold on
errorbar([mean(number_of_registered_neurons(mouse_group==1)),mean(number_of_registered_neurons(mouse_group==2))],[std(number_of_registered_neurons(mouse_group==1))./sqrt(sum(mouse_group==1)),std(number_of_registered_neurons(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',3,'color','k')
hold on
scatter(x_vec_1,number_of_registered_neurons(mouse_group==1),25,'k','filled');
hold on
scatter(x_vec_2,number_of_registered_neurons(mouse_group==2),25,'k','filled');
xlim([0.5 2.5])
ylim([0 1500])
set(gca,'ytick',0:300:1500)
ylabel('Number of registered neurons')
box off
hold on
axis square
set(gca,'xtick',1:3)
set(gca,'xticklabels',{'CA1','CA3'})
set(gca,'fontsize',16)
set(gcf,'PaperPositionMode','auto')
axis square
box off
title('Figure S3D')

% Figure S3E - Percentage of active cells per day from all registered cells:
figure
x_vec_1=1+0.3*(0.5-rand(1,size(fraction_of_neurons(mouse_group==1),2)));
x_vec_2=2+0.3*(0.5-rand(1,size(fraction_of_neurons(mouse_group==2),2)));
bar([mean(fraction_of_neurons(mouse_group==1)),mean(fraction_of_neurons(mouse_group==2))],0.5,'FaceColor','none')
hold on
errorbar([mean(fraction_of_neurons(mouse_group==1)),mean(fraction_of_neurons(mouse_group==2))],[std(fraction_of_neurons(mouse_group==1))./sqrt(sum(mouse_group==1)),std(fraction_of_neurons(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',3,'color','k')
hold on
scatter(x_vec_1,fraction_of_neurons(mouse_group==1),25,'k','filled');
hold on
scatter(x_vec_2,fraction_of_neurons(mouse_group==2),25,'k','filled');
xlim([0.5 2.5])
ylim([0 80])
set(gca,'ytick',0:20:80)
ylabel('% active cells')
box off
hold on
axis square
set(gca,'xtick',1:3)
set(gca,'xticklabels',{'CA1','CA3'})
set(gca,'fontsize',16)
set(gcf,'PaperPositionMode','auto')
axis square
box off
title('Figure S3E')

% Figure S3F - Percentage of overall registration errors (average of false positives and false negatives) per mouse:
figure
x_vec_1=1+0.3*(0.5-rand(1,size(overall_registration_errors_CA1,2)));
x_vec_2=2+0.3*(0.5-rand(1,size(overall_registration_errors_CA3,2)));
bar([mean(overall_registration_errors_CA1),mean(overall_registration_errors_CA3)],0.5,'FaceColor','none')
hold on
errorbar([mean(overall_registration_errors_CA1),mean(overall_registration_errors_CA3)],[std(overall_registration_errors_CA1)./sqrt(length(overall_registration_errors_CA1)),std(overall_registration_errors_CA3)./sqrt(length(overall_registration_errors_CA3))],'.','linewidth',3,'color','k')
hold on
scatter(x_vec_1,overall_registration_errors_CA1,25,'k','filled');
hold on
scatter(x_vec_2,overall_registration_errors_CA3,25,'k','filled');
xlim([0.5 2.5])
ylim([0 30])
set(gca,'ytick',0:10:50)
ylabel('% Registration errors')
box off
hold on
axis square
set(gca,'xtick',1:3)
set(gca,'xticklabels',{'CA1','CA3'})
set(gca,'fontsize',16)
set(gcf,'PaperPositionMode','auto')
axis square
box off
title('Figure S3F')

% Figure S3G-I - Recapitulating the results for a fixed number of registered cells:
% Figure S3G - Percentage of active cells over days from all registered cells:
figure
errorbar(days_vec(1:maximal_number_of_sessions/2),100*mean_fraction_of_neurons_dynamics_CA1(1:maximal_number_of_sessions/2),100*std_fraction_of_neurons_dynamics_CA1(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
hold on
errorbar(days_vec(1:maximal_number_of_sessions/2),100*mean_fraction_of_neurons_dynamics_CA3(1:maximal_number_of_sessions/2),100*std_fraction_of_neurons_dynamics_CA3(1:maximal_number_of_sessions/2)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
errorbar(days_vec(maximal_number_of_sessions/2+1:end),100*mean_fraction_of_neurons_dynamics_CA1(maximal_number_of_sessions/2+1:end),100*std_fraction_of_neurons_dynamics_CA1(maximal_number_of_sessions/2+1:end)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
errorbar(days_vec(maximal_number_of_sessions/2+1:end),100*mean_fraction_of_neurons_dynamics_CA3(maximal_number_of_sessions/2+1:end),100*std_fraction_of_neurons_dynamics_CA3(maximal_number_of_sessions/2+1:end)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
plot([17 17],[0 100],'--','color','k','linewidth',2)
xlim([1 33])
ylim([0 80])
xlabel('Time (days)')
ylabel('% active cells')
set(gca,'fontsize',16)
legend('CA1','CA3')
legend('boxoff')
box off
axis square
title('Figure S3G')


% Figure S3H - Rate correlation versus elapsed time:
figure
errorbar(elapsed_days_vec(2:end),mean_rate_dynamics_CA1(2:end),std_rate_dynamics_CA1(2:end)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
hold on
errorbar(elapsed_days_vec(2:end),mean_rate_dynamics_CA3(2:end),std_rate_dynamics_CA3(2:end)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
errorbar(elapsed_days_vec(1),mean_rate_dynamics_CA1(1),std_rate_dynamics_CA1(1)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
errorbar(elapsed_days_vec(1),mean_rate_dynamics_CA3(1),std_rate_dynamics_CA3(1)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
xlim([-0.5 14])
ylim([0 1])
xlabel('Elapsed time (days)')
ylabel('Rate correlation')
set(gca,'fontsize',16)
legend('CA1','CA3')
legend('boxoff')
box off
axis square
title('Figure S3H - Across days')

% Figure S3I - PV correlation versus elapsed time:
figure
errorbar(elapsed_days_vec(2:end),mean_PV_dynamics_CA1(2:end),std_PV_dynamics_CA1(2:end)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
hold on
errorbar(elapsed_days_vec(2:end),mean_PV_dynamics_CA3(2:end),std_PV_dynamics_CA3(2:end)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
errorbar(elapsed_days_vec(1),mean_PV_dynamics_CA1(1),std_PV_dynamics_CA1(1)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
errorbar(elapsed_days_vec(1),mean_PV_dynamics_CA3(1),std_PV_dynamics_CA3(1)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
xlim([-0.5 14])
ylim([0 1])
set(gca,'ytick',-0.2:0.2:1)
xlabel('Elapsed time (days)')
ylabel('PV correlation')
set(gca,'fontsize',16)
box off
axis square
legend('CA1','CA3')
legend('boxoff')
title('Figure S3I - Across days')

%% Figure S4 - Additional validations of the observed assembly organization in CA3

% Figure S4A-C - Relationship between tuning precision and pairwise correlations and comparison across the two running
% directions:
within_day_pairwise_tuning_corr_single_cells=cell(1,number_of_mice);
SI_SSR_single_cells=cell(1,number_of_mice);

for n=1:number_of_mice
    for k=1:maximal_number_of_sessions
        within_day_pairwise_tuning_corr_single_cells{n}=[within_day_pairwise_tuning_corr_single_cells{n},across_mice_data{n}.within_day_pairwise_tuning_corr_single_cells{k}];
        SI_SSR_single_cells{n}=[SI_SSR_single_cells{n},across_mice_data{n}.within_day_SI_SSR_single_cells{k}];
    end
end

for n=1:number_of_mice
    within_day_pairwise_tuning_corr_single_cells{n}(isnan(SI_SSR_single_cells{n}))=[];
    SI_SSR_single_cells{n}(isnan(SI_SSR_single_cells{n}))=[];
    
    SI_SSR_single_cells{n}(isnan(within_day_pairwise_tuning_corr_single_cells{n}))=[];
    within_day_pairwise_tuning_corr_single_cells{n}(isnan(within_day_pairwise_tuning_corr_single_cells{n}))=[];
end

information_vec=0.4:0.2:3;
average_info_vec=nan(1,length(information_vec)-1);
average_within_pairwise_correlation_per_info=nan(number_of_mice,length(information_vec)-1);
for k=1:number_of_mice
    for n=1:length(information_vec)-1
        this_info_indexes=find(SI_SSR_single_cells{k}>=information_vec(n) & SI_SSR_single_cells{k}<=information_vec(n+1));
        average_within_pairwise_correlation_per_info(k,n)=mean(within_day_pairwise_tuning_corr_single_cells{k}(this_info_indexes),'omitnan');
        if k==1
            average_info_vec(n)=mean(information_vec(n:n+1));
        end
    end
end

mean_within_pairwise_correlation_per_info_CA1=mean(average_within_pairwise_correlation_per_info(mouse_group==1,:),'omitnan');
sd_within_pairwise_correlation_per_info_CA1=std(average_within_pairwise_correlation_per_info(mouse_group==1,:),'omitnan');
mean_within_pairwise_correlation_per_info_CA3=mean(average_within_pairwise_correlation_per_info(mouse_group==2,:),'omitnan');
sd_within_pairwise_correlation_per_info_CA3=std(average_within_pairwise_correlation_per_info(mouse_group==2,:),'omitnan');

mean_SI_SSR_single_cells_separated_novel=nan(number_of_mice,2,length(novel_sessions));
mean_pairwise_corr_single_cells_separated_novel=nan(number_of_mice,2,length(novel_sessions));
for k=1:number_of_mice
    for n=1:length(novel_sessions)
        mean_SI_SSR_single_cells_separated_novel(k,1,n)=mean(across_mice_data{k}.SI_SSR_single_cells_separated{1,novel_sessions(n)},'omitnan');
        mean_SI_SSR_single_cells_separated_novel(k,2,n)=mean(across_mice_data{k}.SI_SSR_single_cells_separated{2,novel_sessions(n)},'omitnan');
        mean_pairwise_corr_single_cells_separated_novel(k,1,n)=mean(across_mice_data{k}.within_day_pairwise_tuning_corr_single_cells_separated{1,novel_sessions(n)},'omitnan');
        mean_pairwise_corr_single_cells_separated_novel(k,2,n)=mean(across_mice_data{k}.within_day_pairwise_tuning_corr_single_cells_separated{2,novel_sessions(n)},'omitnan');
    end
end

% Figure S4D-F - recapitulating the results by comparing against the across
% mice distribution and removing crosstalk:
corr_vec=-0.975:0.05:0.975;

same_bin_across_mice_pairwise_correlations_CA1_novel_AM=[];
same_bin_across_mice_pairwise_correlations_CA3_novel_AM=[];
normalized_same_bin_within_day_pair_corr_per_mouse_novel_AM=nan(number_of_mice,length(corr_vec));
normalized_same_bin_across_mice_pair_corr_per_mouse_novel_AM=nan(number_of_mice,length(corr_vec));
for n=1:number_of_mice
    temp_hist=hist(across_mice_data{n}.same_bin_within_day_pairwise_correlations_per_mouse_novel_AM,corr_vec);
    normalized_temp_hist=temp_hist./sum(temp_hist);
    normalized_same_bin_within_day_pair_corr_per_mouse_novel_AM(n,:)=normalized_temp_hist;
    
    temp_hist=hist(across_mice_data{n}.same_bin_across_mice_pairwise_correlations_per_mouse_novel_AM,corr_vec);
    normalized_temp_hist=temp_hist./sum(temp_hist);
    normalized_same_bin_across_mice_pair_corr_per_mouse_novel_AM(n,:)=normalized_temp_hist;
    if mouse_group(n)==1
        same_bin_across_mice_pairwise_correlations_CA1_novel_AM=[same_bin_across_mice_pairwise_correlations_CA1_novel_AM,across_mice_data{n}.same_bin_across_mice_pairwise_correlations_per_mouse_novel_AM];
    else
        same_bin_across_mice_pairwise_correlations_CA3_novel_AM=[same_bin_across_mice_pairwise_correlations_CA3_novel_AM,across_mice_data{n}.same_bin_across_mice_pairwise_correlations_per_mouse_novel_AM];
    end
end
same_bin_across_mice_pairwise_dist_CA1_novel_AM=hist(same_bin_across_mice_pairwise_correlations_CA1_novel_AM,corr_vec);
normalized_same_bin_across_mice_pairwise_dist_CA1_novel_AM=same_bin_across_mice_pairwise_dist_CA1_novel_AM./sum(same_bin_across_mice_pairwise_dist_CA1_novel_AM);
same_bin_across_mice_pairwise_dist_CA3_novel_AM=hist(same_bin_across_mice_pairwise_correlations_CA3_novel_AM,corr_vec);
normalized_same_bin_across_mice_pairwise_dist_CA3_novel_AM=same_bin_across_mice_pairwise_dist_CA3_novel_AM./sum(same_bin_across_mice_pairwise_dist_CA3_novel_AM);
cdf_same_bin_within_day_pair_corr_per_mouse_novel_AM=cumsum(normalized_same_bin_within_day_pair_corr_per_mouse_novel_AM')';
cdf_same_bin_across_mice_pair_corr_per_mouse_novel_AM=cumsum(normalized_same_bin_across_mice_pair_corr_per_mouse_novel_AM')';

normalized_same_bin_within_day_pair_corr_per_mouse_novel=nan(number_of_mice,length(corr_vec));
normalized_same_bin_across_mice_pair_corr_per_mouse_novel=nan(number_of_mice,length(corr_vec));
normalized_same_bin_across_sides_pair_corr_per_mouse_novel=nan(number_of_mice,length(corr_vec));
for n=1:number_of_mice
    temp_hist=hist(across_mice_data{n}.same_bin_within_day_pairwise_correlations_per_mouse_novel_CTC,corr_vec);
    normalized_temp_hist=temp_hist./sum(temp_hist);
    normalized_same_bin_within_day_pair_corr_per_mouse_novel(n,:)=normalized_temp_hist;
    
    temp_hist=hist(across_mice_data{n}.same_bin_across_mice_pairwise_correlations_per_mouse_novel_CTC,corr_vec);
    normalized_temp_hist=temp_hist./sum(temp_hist);
    normalized_same_bin_across_mice_pair_corr_per_mouse_novel(n,:)=normalized_temp_hist;
    
    temp_hist=hist(across_mice_data{n}.same_bin_across_sides_pairwise_correlations_per_mouse_novel_CTC,corr_vec);
    normalized_temp_hist=temp_hist./sum(temp_hist);
    normalized_same_bin_across_sides_pair_corr_per_mouse_novel(n,:)=normalized_temp_hist;
end

cdf_same_bin_within_day_pair_corr_per_mouse_novel=cumsum(normalized_same_bin_within_day_pair_corr_per_mouse_novel')';
cdf_same_bin_across_mice_pair_corr_per_mouse_novel=cumsum(normalized_same_bin_across_mice_pair_corr_per_mouse_novel')';
cdf_same_bin_across_sides_pair_corr_per_mouse_novel=cumsum(normalized_same_bin_across_sides_pair_corr_per_mouse_novel')';

% Figure S4G-I - Maintenance of tuning properties across environments:
maximal_PV_matrix=nan(number_of_mice,2*maximal_number_of_sessions,2*maximal_number_of_sessions);

for n=1:number_of_mice
    maximal_PV_matrix(n,:,:)=across_mice_data{n}.maximal_PV_correlations;
end

for n=1:number_of_mice
    for k=1:2*maximal_number_of_sessions
        maximal_PV_matrix(n,k,k)=1;
    end
end
across_env_maximal_PV_correlations=mean(squeeze(mean(mean(maximal_PV_matrix(:,15:16,17:18),2,'omitnan'),2,'omitnan')),2,'omitnan');

all_pairwise_correlations_session_one_right_across=cell(1,number_of_mice);
all_pairwise_correlations_session_one_left_across=cell(1,number_of_mice);
all_pairwise_correlations_session_one_right_within=cell(1,number_of_mice);
all_pairwise_correlations_session_one_left_within=cell(1,number_of_mice);
all_pairwise_correlations_session_two_right_across=cell(1,number_of_mice);
all_pairwise_correlations_session_two_left_across=cell(1,number_of_mice);
all_pairwise_correlations_session_two_right_within=cell(1,number_of_mice);
all_pairwise_correlations_session_two_left_within=cell(1,number_of_mice);
for n=1:number_of_mice
    all_pairwise_correlations_session_one_right_across{n}=across_mice_data{n}.all_pairwise_correlations_session_one_right_across;
    all_pairwise_correlations_session_one_left_across{n}=across_mice_data{n}.all_pairwise_correlations_session_one_left_across;
    all_pairwise_correlations_session_one_right_within{n}=across_mice_data{n}.all_pairwise_correlations_session_one_right_within;
    all_pairwise_correlations_session_one_left_within{n}=across_mice_data{n}.all_pairwise_correlations_session_one_left_within;
    all_pairwise_correlations_session_two_right_across{n}=across_mice_data{n}.all_pairwise_correlations_session_two_right_across;
    all_pairwise_correlations_session_two_left_across{n}=across_mice_data{n}.all_pairwise_correlations_session_two_left_across;
    all_pairwise_correlations_session_two_right_within{n}=across_mice_data{n}.all_pairwise_correlations_session_two_right_within;
    all_pairwise_correlations_session_two_left_within{n}=across_mice_data{n}.all_pairwise_correlations_session_two_left_within;
end

pairwise_correlation_vec=-0.5:0.15:1;
average_pairwise_correlation_vec=nan(1,length(pairwise_correlation_vec)-1);
average_pairwise_correlation_envA_per_correlation_envA=nan(number_of_mice,length(pairwise_correlation_vec)-1);
average_pairwise_correlation_envA_per_correlation_envB=nan(number_of_mice,length(pairwise_correlation_vec)-1);
for k=1:number_of_mice
    for n=1:length(pairwise_correlation_vec)-1
        this_pairwise_correlation_indexes=find(all_pairwise_correlations_session_one_right_within{k}>=pairwise_correlation_vec(n) & all_pairwise_correlations_session_one_right_within{k}<=pairwise_correlation_vec(n+1));
        temp_correlations_right=all_pairwise_correlations_session_two_right_within{k}(this_pairwise_correlation_indexes);
        this_pairwise_correlation_indexes=find(all_pairwise_correlations_session_one_left_within{k}>=pairwise_correlation_vec(n) & all_pairwise_correlations_session_one_left_within{k}<=pairwise_correlation_vec(n+1));
        temp_correlations_left=all_pairwise_correlations_session_two_left_within{k}(this_pairwise_correlation_indexes);
        average_pairwise_correlation_envA_per_correlation_envA(k,n)=mean([temp_correlations_right,temp_correlations_left],'omitnan');
        this_pairwise_correlation_indexes= all_pairwise_correlations_session_one_right_across{k}>=pairwise_correlation_vec(n) & all_pairwise_correlations_session_one_right_across{k}<=pairwise_correlation_vec(n+1);
        temp_correlations_right=all_pairwise_correlations_session_two_right_across{k}(this_pairwise_correlation_indexes);
        this_pairwise_correlation_indexes=find(all_pairwise_correlations_session_one_left_across{k}>=pairwise_correlation_vec(n) & all_pairwise_correlations_session_one_left_across{k}<=pairwise_correlation_vec(n+1));
        temp_correlations_left=all_pairwise_correlations_session_two_left_across{k}(this_pairwise_correlation_indexes);
        average_pairwise_correlation_envA_per_correlation_envB(k,n)=mean([temp_correlations_right,temp_correlations_left],'omitnan');
        if k==1
            average_pairwise_correlation_vec(n)=mean(pairwise_correlation_vec(n:n+1));
        end
    end
end

mean_pairwise_correlation_vec_within_CA1=mean(average_pairwise_correlation_envA_per_correlation_envA(mouse_group==1,:),'omitnan');
sd_pairwise_correlation_vec_within_CA1=std(average_pairwise_correlation_envA_per_correlation_envA(mouse_group==1,:),'omitnan');
mean_pairwise_correlation_vec_within_CA3=mean(average_pairwise_correlation_envA_per_correlation_envA(mouse_group==2,:),'omitnan');
sd_pairwise_correlation_vec_within_CA3=std(average_pairwise_correlation_envA_per_correlation_envA(mouse_group==2,:),'omitnan');
mean_pairwise_correlation_vec_across_CA1=mean(average_pairwise_correlation_envA_per_correlation_envB(mouse_group==1,:),'omitnan');
sd_pairwise_correlation_vec_across_CA1=std(average_pairwise_correlation_envA_per_correlation_envB(mouse_group==1,:),'omitnan');
mean_pairwise_correlation_vec_across_CA3=mean(average_pairwise_correlation_envA_per_correlation_envB(mouse_group==2,:),'omitnan');
sd_pairwise_correlation_vec_across_CA3=std(average_pairwise_correlation_envA_per_correlation_envB(mouse_group==2,:),'omitnan');

significant_cells_is_sig_dynamics=nan(number_of_mice,maximal_number_of_sessions);
not_significant_cells_is_sig_dynamics=nan(number_of_mice,maximal_number_of_sessions);

significant_cells_is_sig_across_env=nan(number_of_mice,maximal_number_of_sessions);
not_significant_cells_is_sig_across_env=nan(number_of_mice,maximal_number_of_sessions);

for n=1:number_of_mice
    significant_cells_is_sig_dynamics(n,:)=across_mice_data{n}.mean_significant_cells_is_sig_dynamics;
    not_significant_cells_is_sig_dynamics(n,:)=across_mice_data{n}.mean_not_significant_cells_is_sig_dynamics;
    significant_cells_is_sig_across_env(n,:)=across_mice_data{n}.mean_significant_cells_is_sig_across_env;
    not_significant_cells_is_sig_across_env(n,:)=across_mice_data{n}.mean_not_significant_cells_is_sig_across_env;
end

mean_significant_cells_is_sig_across_env=significant_cells_is_sig_across_env(:,2);
mean_not_significant_cells_is_sig_across_env=not_significant_cells_is_sig_across_env(:,2);
mean_significant_cells_is_sig_within_env=significant_cells_is_sig_dynamics(:,2);
mean_not_significant_cells_is_sig_within_env=not_significant_cells_is_sig_dynamics(:,2);

% Figure S4J-L - Analysis of pairwise noise correlations:
noise_correlations_same_pos_correlated_single_cells_dynamics=cell(number_of_mice,maximal_number_of_sessions);
noise_correlations_same_pos_uncorrelated_single_cells_dynamics=cell(number_of_mice,maximal_number_of_sessions);
for n=1:number_of_mice
    for k=1:maximal_number_of_sessions
        noise_correlations_same_pos_correlated_single_cells_dynamics{n,k}=across_mice_data{n}.noise_correlations_same_pos_correlated_single_cells{k};
        noise_correlations_same_pos_uncorrelated_single_cells_dynamics{n,k}=across_mice_data{n}.noise_correlations_same_pos_uncorrelated_single_cells{k};
    end
end

same_pos_noise_correlations_per_mouse_novel=cell(1,number_of_mice);
for n=1:number_of_mice
    for k=1:maximal_number_of_sessions
        if ~isempty(intersect(k,novel_sessions))
            same_pos_noise_correlations_per_mouse_novel{n}=[same_pos_noise_correlations_per_mouse_novel{n},noise_correlations_same_pos_correlated_single_cells_dynamics{n,k}];
            same_pos_noise_correlations_per_mouse_novel{n}=[same_pos_noise_correlations_per_mouse_novel{n},noise_correlations_same_pos_uncorrelated_single_cells_dynamics{n,k}];
        end
    end
end

mean_same_pos_noise_correlations_per_mouse_novel=nan(1,number_of_mice);
for n=1:number_of_mice
    mean_same_pos_noise_correlations_per_mouse_novel(n)=mean(same_pos_noise_correlations_per_mouse_novel{n},'omitnan');
end

noise_corr_versus_shift_novel=nan(number_of_mice,num_bins);
shuffle_noise_corr_versus_shift_novel=nan(number_of_mice,num_bins);
for n=1:number_of_mice
    noise_corr_versus_shift_novel(n,:)=across_mice_data{n}.average_noise_corr_versus_shift_novel;
    for k=1:num_bins
        shuffle_noise_corr_versus_shift_novel(n,k)=mean(across_mice_data{n}.shuffle_noise_correlation_versus_shift_novel{k});
    end
end

mean_noise_corr_versus_shift_CA1_novel=mean(noise_corr_versus_shift_novel(mouse_group==1,:),'omitnan');
std_noise_corr_versus_shift_CA1_novel=std(noise_corr_versus_shift_novel(mouse_group==1,:),'omitnan');
mean_noise_corr_versus_shift_CA3_novel=mean(noise_corr_versus_shift_novel(mouse_group==2,:),'omitnan');
std_noise_corr_versus_shift_CA3_novel=std(noise_corr_versus_shift_novel(mouse_group==2,:),'omitnan');

mean_shuffle_noise_corr_versus_shift_novel=mean(shuffle_noise_corr_versus_shift_novel(mouse_group==1 | mouse_group==2,:),'omitnan');
std_shuffle_noise_corr_versus_shift_novel=std(shuffle_noise_corr_versus_shift_novel(mouse_group==1 | mouse_group==2,:),'omitnan');

distance_vector=0:20:500;
average_noise_correlation_versus_distance_novel=nan(number_of_mice,length(distance_vector)-1);
average_shuffle_noise_correlation_versus_distance_novel=nan(number_of_mice,length(distance_vector)-1);

for n=1:number_of_mice
    average_noise_correlation_versus_distance_novel(n,:)=across_mice_data{n}.average_noise_correlation_versus_distance_novel;
    average_shuffle_noise_correlation_versus_distance_novel(n,:)=across_mice_data{n}.average_shuffle_noise_correlation_versus_distance_novel;
end

mean_noise_correlation_versus_distance_CA1_novel=mean(average_noise_correlation_versus_distance_novel(mouse_group==1,:),'omitnan');
std_noise_correlation_versus_distance_CA1_novel=std(average_noise_correlation_versus_distance_novel(mouse_group==1,:),'omitnan');
mean_noise_correlation_versus_distance_CA3_novel=mean(average_noise_correlation_versus_distance_novel(mouse_group==2,:),'omitnan');
std_noise_correlation_versus_distance_CA3_novel=std(average_noise_correlation_versus_distance_novel(mouse_group==2,:),'omitnan');

mean_shuffle_noise_correlation_versus_distance_novel=mean(average_shuffle_noise_correlation_versus_distance_novel,'omitnan');
std_shuffle_noise_correlation_versus_distance_novel=std(average_shuffle_noise_correlation_versus_distance_novel,'omitnan');

% Figure S4A - Pairwise tuning curve correlation versus spatial
% information:
figure
errorbar(average_info_vec,mean_within_pairwise_correlation_per_info_CA1,sd_within_pairwise_correlation_per_info_CA1./sqrt(sum(mouse_group==1)),'-*','markersize',5,'color','b','linewidth',2)
hold on
errorbar(average_info_vec,mean_within_pairwise_correlation_per_info_CA3,sd_within_pairwise_correlation_per_info_CA3./sqrt(sum(mouse_group==2)),'-*','markersize',5,'color','r','linewidth',2)
xlim([0.4 3])
ylim([0.5 1])
xlabel('Spatial information (bit/spike)')
ylabel('Pairwise tuning correlation')
box off
axis square
set(gca,'fontsize',16)
title('Figure S4A')
legend({'CA1','CA3'},'Location','Northwest')
legend('boxoff')

% Figure S4B - Average spatial information compared between the two running directions:
lgd = [];
figure
plot([0 2],[0 2],'--k','linewidth',2)
hold on
for n=1:number_of_mice
    if mouse_group(n)==1
       lgd(1) =  plot(mean(mean_SI_SSR_single_cells_separated_novel(n,1,:),'omitnan'),mean(mean_SI_SSR_single_cells_separated_novel(n,2,:),'omitnan'),'.','markersize',20,'color','b');
        errorbar(mean(mean_SI_SSR_single_cells_separated_novel(n,1,:),'omitnan'),mean(mean_SI_SSR_single_cells_separated_novel(n,2,:),'omitnan'),std(mean_SI_SSR_single_cells_separated_novel(n,2,:),'omitnan'),'color','b')
        errorbar(mean(mean_SI_SSR_single_cells_separated_novel(n,1,:),'omitnan'),mean(mean_SI_SSR_single_cells_separated_novel(n,2,:),'omitnan'),std(mean_SI_SSR_single_cells_separated_novel(n,2,:),'omitnan'),'horizontal','color','b')
    elseif mouse_group(n)==2
        lgd(2) =  plot(mean(mean_SI_SSR_single_cells_separated_novel(n,1,:),'omitnan'),mean(mean_SI_SSR_single_cells_separated_novel(n,2,:),'omitnan'),'.','markersize',20,'color','r');
        errorbar(mean(mean_SI_SSR_single_cells_separated_novel(n,1,:),'omitnan'),mean(mean_SI_SSR_single_cells_separated_novel(n,2,:),'omitnan'),std(mean_SI_SSR_single_cells_separated_novel(n,2,:),'omitnan'),'color','r')
        errorbar(mean(mean_SI_SSR_single_cells_separated_novel(n,1,:),'omitnan'),mean(mean_SI_SSR_single_cells_separated_novel(n,2,:),'omitnan'),std(mean_SI_SSR_single_cells_separated_novel(n,2,:),'omitnan'),'horizontal','color','r')
    end
end
xlabel('Spatial information - right (bit/spike)')
ylabel('Spatial information - left (bit/spike)')
xlim([0 2])
ylim([0 2])
set(gca,'fontsize',16)
legend(lgd,'CA1','CA3','Location','Northwest')
legend boxoff
axis square
box off
title('Figure S4B')

% Figure S4C - Average pairwise correlation compared between the two running directions:
figure
plot([0 2],[0 2],'--k','linewidth',2)
hold on
for n=1:number_of_mice
    if mouse_group(n)==1
        lgd(1) = plot(mean(mean_pairwise_corr_single_cells_separated_novel(n,1,:),'omitnan'),mean(mean_pairwise_corr_single_cells_separated_novel(n,2,:),'omitnan'),'.','markersize',20,'color','b');
        errorbar(mean(mean_pairwise_corr_single_cells_separated_novel(n,1,:),'omitnan'),mean(mean_pairwise_corr_single_cells_separated_novel(n,2,:),'omitnan'),std(mean_pairwise_corr_single_cells_separated_novel(n,2,:),'omitnan'),'color','b')
        errorbar(mean(mean_pairwise_corr_single_cells_separated_novel(n,1,:),'omitnan'),mean(mean_pairwise_corr_single_cells_separated_novel(n,2,:),'omitnan'),std(mean_pairwise_corr_single_cells_separated_novel(n,2,:),'omitnan'),'horizontal','color','b')
    elseif mouse_group(n)==2
        lgd(2) = plot(mean(mean_pairwise_corr_single_cells_separated_novel(n,1,:),'omitnan'),mean(mean_pairwise_corr_single_cells_separated_novel(n,2,:),'omitnan'),'.','markersize',20,'color','r');
        errorbar(mean(mean_pairwise_corr_single_cells_separated_novel(n,1,:),'omitnan'),mean(mean_pairwise_corr_single_cells_separated_novel(n,2,:),'omitnan'),std(mean_pairwise_corr_single_cells_separated_novel(n,2,:),'omitnan'),'color','r')
        errorbar(mean(mean_pairwise_corr_single_cells_separated_novel(n,1,:),'omitnan'),mean(mean_pairwise_corr_single_cells_separated_novel(n,2,:),'omitnan'),std(mean_pairwise_corr_single_cells_separated_novel(n,2,:),'omitnan'),'horizontal','color','r')
    end
end
xlabel('Pairwise tuning correlation - right')
ylabel('Pairwise tuning correlation - left')
xlim([0.5 1])
ylim([0.5 1])
set(gca,'fontsize',16)
axis square
box off
legend(lgd,'CA1','CA3','Location','Northwest')
legend boxoff
title('Figure S4C')

% Figure S4D - Comparing to a null distribution obtained from different mice:
mean_cdf_same_bin_within_day_pair_corr_CA1_novel=mean(cdf_same_bin_within_day_pair_corr_per_mouse_novel_AM(mouse_group==1,:));
sd_cdf_same_bin_within_day_pair_corr_CA1_novel=std(cdf_same_bin_within_day_pair_corr_per_mouse_novel_AM(mouse_group==1,:));
mean_cdf_same_bin_within_day_pair_corr_CA3_novel=mean(cdf_same_bin_within_day_pair_corr_per_mouse_novel_AM(mouse_group==2,:));
sd_cdf_same_bin_within_day_pair_corr_CA3_novel=std(cdf_same_bin_within_day_pair_corr_per_mouse_novel_AM(mouse_group==2,:));

figure
errorbar(corr_vec,1-fliplr(mean_cdf_same_bin_within_day_pair_corr_CA1_novel),fliplr(sd_cdf_same_bin_within_day_pair_corr_CA1_novel)./sqrt(sum(mouse_group==1)),'-ob','markersize',5,'linewidth',1)
hold on
plot(corr_vec,1-fliplr(cumsum(normalized_same_bin_across_mice_pairwise_dist_CA1_novel_AM)),'--b','linewidth',2)
errorbar(corr_vec,1-fliplr(mean_cdf_same_bin_within_day_pair_corr_CA3_novel),fliplr(sd_cdf_same_bin_within_day_pair_corr_CA3_novel)./sqrt(sum(mouse_group==1)),'-or','markersize',5,'linewidth',1)
plot(corr_vec,1-fliplr(cumsum(normalized_same_bin_across_mice_pairwise_dist_CA3_novel_AM)),'--r','linewidth',2)
errorbar(corr_vec,1-fliplr(mean_cdf_same_bin_within_day_pair_corr_CA1_novel),fliplr(sd_cdf_same_bin_within_day_pair_corr_CA1_novel)./sqrt(sum(mouse_group==1)),'-ob','markersize',5,'linewidth',1)
xlim([-1 0])
ylim([0 1])
xlabel('Tuning curve correlation')
ylabel('Cumulative fraction of cell-pairs')
set(gca,'xtick',-1:0.5:1)
set(gca,'xticklabels',1:-0.5:-1)
legend('CA1 - Within day','CA1 - Across mice','CA3 - Within day','CA3 - Across mice','Location','Southeast')
legend boxoff
axis square
set(gca,'fontsize',16)
box off
title('Figure S4D - Within vs. across mice')

% Figure S4E - Comparing to a null distribution obtained from different mice:
cdf_diff_same_bin_within_day_pair_corr_CA1_novel=(1-fliplr(cdf_same_bin_within_day_pair_corr_per_mouse_novel_AM(mouse_group==1,:)))-(1-fliplr(cdf_same_bin_across_mice_pair_corr_per_mouse_novel(mouse_group==1,:)));
mean_cdf_diff_same_bin_within_day_pair_corr_CA1_novel=mean(cdf_diff_same_bin_within_day_pair_corr_CA1_novel);
sd_cdf_diff_same_bin_within_day_pair_corr_CA1_novel=std(cdf_diff_same_bin_within_day_pair_corr_CA1_novel);

cdf_diff_same_bin_within_day_pair_corr_CA3_novel=(1-fliplr(cdf_same_bin_within_day_pair_corr_per_mouse_novel_AM(mouse_group==2,:)))-(1-fliplr(cdf_same_bin_across_mice_pair_corr_per_mouse_novel(mouse_group==2,:)));
mean_cdf_diff_same_bin_within_day_pair_corr_CA3_novel=mean(cdf_diff_same_bin_within_day_pair_corr_CA3_novel);
sd_cdf_diff_same_bin_within_day_pair_corr_CA3_novel=std(cdf_diff_same_bin_within_day_pair_corr_CA3_novel);

figure
errorbar(corr_vec,mean_cdf_diff_same_bin_within_day_pair_corr_CA1_novel,sd_cdf_diff_same_bin_within_day_pair_corr_CA1_novel./sqrt(sum(mouse_group==1)),'-b','linewidth',2)
hold on
errorbar(corr_vec,mean_cdf_diff_same_bin_within_day_pair_corr_CA3_novel,sd_cdf_diff_same_bin_within_day_pair_corr_CA3_novel./sqrt(sum(mouse_group==2)),'-r','linewidth',2)
plot([-1 0],[0 0],'--k','linewidth',2)
xlim([-1 0])
ylim([-0.05 0.15])
xlabel('Tuning curve correlation')
ylabel('Cumulative fraction difference')
set(gca,'xtick',-1:0.5:1)
set(gca,'xticklabels',1:-0.5:-1)
legend('CA1','CA3')
legend boxoff
axis square
set(gca,'fontsize',16)
box off
title('Figure S4E - Within vs. across mice')

% Figure S4F - Removing crosstalk by excluding cells that are < 50 microns
% apart (comparison with across mice):
cdf_diff_same_bin_within_day_pair_corr_CA1_novel=(1-fliplr(cdf_same_bin_within_day_pair_corr_per_mouse_novel(mouse_group==1,:)))-(1-fliplr(cdf_same_bin_across_mice_pair_corr_per_mouse_novel(mouse_group==1,:)));
mean_cdf_diff_same_bin_within_day_pair_corr_CA1_novel=mean(cdf_diff_same_bin_within_day_pair_corr_CA1_novel);
sd_cdf_diff_same_bin_within_day_pair_corr_CA1_novel=std(cdf_diff_same_bin_within_day_pair_corr_CA1_novel);

cdf_diff_same_bin_within_day_pair_corr_CA3_novel=(1-fliplr(cdf_same_bin_within_day_pair_corr_per_mouse_novel(mouse_group==2,:)))-(1-fliplr(cdf_same_bin_across_mice_pair_corr_per_mouse_novel(mouse_group==2,:)));
mean_cdf_diff_same_bin_within_day_pair_corr_CA3_novel=mean(cdf_diff_same_bin_within_day_pair_corr_CA3_novel);
sd_cdf_diff_same_bin_within_day_pair_corr_CA3_novel=std(cdf_diff_same_bin_within_day_pair_corr_CA3_novel);

figure
errorbar(corr_vec,mean_cdf_diff_same_bin_within_day_pair_corr_CA1_novel,sd_cdf_diff_same_bin_within_day_pair_corr_CA1_novel./sqrt(sum(mouse_group==1)),'-b','linewidth',2)
hold on
errorbar(corr_vec,mean_cdf_diff_same_bin_within_day_pair_corr_CA3_novel,sd_cdf_diff_same_bin_within_day_pair_corr_CA3_novel./sqrt(sum(mouse_group==2)),'-r','linewidth',2)
plot([-1 0],[0 0],'--k','linewidth',2)
xlim([-1 0])
ylim([-0.05 0.15])
xlabel('Tuning curve correlation')
ylabel('Cumulative fraction difference')
set(gca,'xtick',-1:0.5:1)
set(gca,'xticklabels',1:-0.5:-1)
legend('CA1','CA3')
legend boxoff
axis square
set(gca,'fontsize',16)
box off
title('Figure S4F - Across mice > 50 microns apart')

% Figure S4F, inset - Removing crosstalk by excluding cells that are < 50
% microns apart (comparison with across sides):
cdf_diff_same_bin_within_day_pair_corr_CA1_novel=(1-fliplr(cdf_same_bin_within_day_pair_corr_per_mouse_novel(mouse_group==1,:)))-(1-fliplr(cdf_same_bin_across_sides_pair_corr_per_mouse_novel(mouse_group==1,:)));
mean_cdf_diff_same_bin_within_day_pair_corr_CA1_novel=mean(cdf_diff_same_bin_within_day_pair_corr_CA1_novel);
sd_cdf_diff_same_bin_within_day_pair_corr_CA1_novel=std(cdf_diff_same_bin_within_day_pair_corr_CA1_novel);

cdf_diff_same_bin_within_day_pair_corr_CA3_novel=(1-fliplr(cdf_same_bin_within_day_pair_corr_per_mouse_novel(mouse_group==2,:)))-(1-fliplr(cdf_same_bin_across_sides_pair_corr_per_mouse_novel(mouse_group==2,:)));
mean_cdf_diff_same_bin_within_day_pair_corr_CA3_novel=mean(cdf_diff_same_bin_within_day_pair_corr_CA3_novel);
sd_cdf_diff_same_bin_within_day_pair_corr_CA3_novel=std(cdf_diff_same_bin_within_day_pair_corr_CA3_novel);

figure
errorbar(corr_vec,mean_cdf_diff_same_bin_within_day_pair_corr_CA1_novel,sd_cdf_diff_same_bin_within_day_pair_corr_CA1_novel./sqrt(sum(mouse_group==1)),'-b','linewidth',2)
hold on
errorbar(corr_vec,mean_cdf_diff_same_bin_within_day_pair_corr_CA3_novel,sd_cdf_diff_same_bin_within_day_pair_corr_CA3_novel./sqrt(sum(mouse_group==2)),'-r','linewidth',2)
plot([-1 0],[0 0],'--k','linewidth',2)
xlim([-1 0])
ylim([-0.05 0.1])
xlabel('Tuning curve correlation')
ylabel('Cumulative fraction difference')
set(gca,'xtick',-1:0.5:1)
set(gca,'xticklabels',1:-0.5:-1)
legend('CA1','CA3')
legend boxoff
axis square
set(gca,'fontsize',16)
box off
title('Figure S4F, inset - Across sides > 50 microns apart')

% Figure S4G - Maintenance of PV correlations across environments:
x_vec_1=1+0.3*(0.5-rand(1,length(across_env_maximal_PV_correlations(mouse_group==1))));
x_vec_2=2+0.3*(0.5-rand(1,length(across_env_maximal_PV_correlations(mouse_group==2))));
figure
bar([1,2],[mean(across_env_maximal_PV_correlations(mouse_group==1)),mean(across_env_maximal_PV_correlations(mouse_group==2))],0.5,'FaceColor','none')
hold on
errorbar([1,2],[mean(across_env_maximal_PV_correlations(mouse_group==1)),mean(across_env_maximal_PV_correlations(mouse_group==2))],[std(across_env_maximal_PV_correlations(mouse_group==1))./sqrt(sum(mouse_group==1)),std(across_env_maximal_PV_correlations(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',3,'color','k')
scatter(x_vec_1,across_env_maximal_PV_correlations(mouse_group==1),15,'k','filled');
scatter(x_vec_2,across_env_maximal_PV_correlations(mouse_group==2),15,'k','filled');
xlim([0.5 2.5])
ylim([-0.1 0.3])
set(gca,'ytick',0:0.1:0.3)
ylabel('PV correlation')
set(gca,'fontsize',16)
box off
hold on
axis square
set(gca,'xtick',[1,2,4,5])
set(gca,'xticklabels',{'CA1','CA3'})
set(gca,'fontsize',16)
set(gcf,'PaperPositionMode','auto')
axis square
box off
title('Figure S4G')

% Figure S4H - Maintenance of pairwise tuning curve correlations across days and environments:
figure
errorbar(average_pairwise_correlation_vec,mean_pairwise_correlation_vec_across_CA1,sd_pairwise_correlation_vec_across_CA1./sqrt(sum(mouse_group==1)),'-b','linewidth',2)
hold on
errorbar(average_pairwise_correlation_vec,mean_pairwise_correlation_vec_across_CA3,sd_pairwise_correlation_vec_across_CA3./sqrt(sum(mouse_group==2)),'-r','linewidth',2)
errorbar(average_pairwise_correlation_vec,mean_pairwise_correlation_vec_within_CA1,sd_pairwise_correlation_vec_within_CA1./sqrt(sum(mouse_group==1)),'--b','linewidth',2)
errorbar(average_pairwise_correlation_vec,mean_pairwise_correlation_vec_within_CA3,sd_pairwise_correlation_vec_within_CA3./sqrt(sum(mouse_group==2)),'--r','linewidth',2)
plot([-0.5 1],[0 0],'--k','linewidth',2)
ylim([-0.5 1])
xlabel('Pairwise correlation - day 1')
ylabel('Pairwise correlation - day 2')
set(gca,'fontsize',16)
legend('CA1 - Across environments','CA3 - Across environments','CA1 - Within environment','CA3 - Within environment')
legend boxoff
axis square
box off
title('Figure S4H')

% Figure S4I -  Across environments recurrence probability for peer-dependent cells:
figure
x_vec_1=1+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_2=2+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_3=4+0.3*(0.5-rand(1,sum(mouse_group==2)));
x_vec_4=5+0.3*(0.5-rand(1,sum(mouse_group==2)));
plot([1 2],[mean(mean_significant_cells_is_sig_across_env(mouse_group==1)), mean(mean_not_significant_cells_is_sig_across_env(mouse_group==1))],'.-b','markersize',30,'linewidth',2)
hold on
plot([4 5],[mean(mean_significant_cells_is_sig_across_env(mouse_group==2)), mean(mean_not_significant_cells_is_sig_across_env(mouse_group==2))],'.-r','markersize',30,'linewidth',2)
errorbar([1 2],[mean(mean_significant_cells_is_sig_across_env(mouse_group==1)),mean(mean_not_significant_cells_is_sig_across_env(mouse_group==1))],[std(mean_significant_cells_is_sig_across_env(mouse_group==1))./sqrt(sum(mouse_group==1)),std(mean_not_significant_cells_is_sig_across_env(mouse_group==1))./sqrt(sum(mouse_group==1))],'.','linewidth',2,'color','b')
errorbar([4 5],[mean(mean_significant_cells_is_sig_across_env(mouse_group==2)),mean(mean_not_significant_cells_is_sig_across_env(mouse_group==2))],[std(mean_significant_cells_is_sig_across_env(mouse_group==2))./sqrt(sum(mouse_group==2)),std(mean_not_significant_cells_is_sig_across_env(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','r')
for n=1:sum(mouse_group==1)
    plot([x_vec_1(n), x_vec_2(n)],[mean_significant_cells_is_sig_across_env(CA1_indexes(n)),mean_not_significant_cells_is_sig_across_env(CA1_indexes(n))],'-','color',[0.7 0.7 0.7])
end
for n=1:sum(mouse_group==2)
    plot([x_vec_3(n), x_vec_4(n)],[mean_significant_cells_is_sig_across_env(CA3_indexes(n)),mean_not_significant_cells_is_sig_across_env(CA3_indexes(n))],'-','color',[0.7 0.7 0.7])
end
plot([1 2],[mean(mean_significant_cells_is_sig_across_env(mouse_group==1)), mean(mean_not_significant_cells_is_sig_across_env(mouse_group==1))],'.-b','markersize',30,'linewidth',2)
plot([4 5],[mean(mean_significant_cells_is_sig_across_env(mouse_group==2)), mean(mean_not_significant_cells_is_sig_across_env(mouse_group==2))],'.-r','markersize',30,'linewidth',2)
errorbar([1 2],[mean(mean_significant_cells_is_sig_across_env(mouse_group==1)),mean(mean_not_significant_cells_is_sig_across_env(mouse_group==1))],[std(mean_significant_cells_is_sig_across_env(mouse_group==1))./sqrt(sum(mouse_group==1)),std(mean_not_significant_cells_is_sig_across_env(mouse_group==1))./sqrt(sum(mouse_group==1))],'.','linewidth',2,'color','b')
errorbar([4 5],[mean(mean_significant_cells_is_sig_across_env(mouse_group==2)),mean(mean_not_significant_cells_is_sig_across_env(mouse_group==2))],[std(mean_significant_cells_is_sig_across_env(mouse_group==2))./sqrt(sum(mouse_group==2)),std(mean_not_significant_cells_is_sig_across_env(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','r')
xlim([0.5 5.5])
ylim([0 0.3])
ylabel('Across env. recurrence probability')
box off
hold on
axis square
set(gca,'xtick',[1 2 4 5])
set(gca,'xticklabels',{'Dep','Ind','Dep','Ind'})
set(gca,'ytick',0:0.1:0.3)
set(gca,'fontsize',16)
axis square
box off
title('Figure S4I')

% Figure S4I, inset -  Within environment recurrence probability for peer-dependent cells:
figure
x_vec_1=1+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_2=2+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_3=4+0.3*(0.5-rand(1,sum(mouse_group==2)));
x_vec_4=5+0.3*(0.5-rand(1,sum(mouse_group==2)));
plot([1 2],[mean(mean_significant_cells_is_sig_within_env(mouse_group==1)), mean(mean_not_significant_cells_is_sig_within_env(mouse_group==1))],'.-b','markersize',30,'linewidth',2)
hold on
plot([4 5],[mean(mean_significant_cells_is_sig_within_env(mouse_group==2)), mean(mean_not_significant_cells_is_sig_within_env(mouse_group==2))],'.-r','markersize',30,'linewidth',2)
errorbar([1 2],[mean(mean_significant_cells_is_sig_within_env(mouse_group==1)),mean(mean_not_significant_cells_is_sig_within_env(mouse_group==1))],[std(mean_significant_cells_is_sig_within_env(mouse_group==1))./sqrt(sum(mouse_group==1)),std(mean_not_significant_cells_is_sig_within_env(mouse_group==1))./sqrt(sum(mouse_group==1))],'.','linewidth',2,'color','b')
errorbar([4 5],[mean(mean_significant_cells_is_sig_within_env(mouse_group==2)),mean(mean_not_significant_cells_is_sig_within_env(mouse_group==2))],[std(mean_significant_cells_is_sig_within_env(mouse_group==2))./sqrt(sum(mouse_group==2)),std(mean_not_significant_cells_is_sig_within_env(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','r')
for n=1:sum(mouse_group==1)
    plot([x_vec_1(n), x_vec_2(n)],[mean_significant_cells_is_sig_within_env(CA1_indexes(n)),mean_not_significant_cells_is_sig_within_env(CA1_indexes(n))],'-','color',[0.7 0.7 0.7])
end
for n=1:sum(mouse_group==2)
    plot([x_vec_3(n), x_vec_4(n)],[mean_significant_cells_is_sig_within_env(CA3_indexes(n)),mean_not_significant_cells_is_sig_within_env(CA3_indexes(n))],'-','color',[0.7 0.7 0.7])
end
plot([1 2],[mean(mean_significant_cells_is_sig_within_env(mouse_group==1)), mean(mean_not_significant_cells_is_sig_within_env(mouse_group==1))],'.-b','markersize',30,'linewidth',2)
plot([4 5],[mean(mean_significant_cells_is_sig_within_env(mouse_group==2)), mean(mean_not_significant_cells_is_sig_within_env(mouse_group==2))],'.-r','markersize',30,'linewidth',2)
errorbar([1 2],[mean(mean_significant_cells_is_sig_within_env(mouse_group==1)),mean(mean_not_significant_cells_is_sig_within_env(mouse_group==1))],[std(mean_significant_cells_is_sig_within_env(mouse_group==1))./sqrt(sum(mouse_group==1)),std(mean_not_significant_cells_is_sig_within_env(mouse_group==1))./sqrt(sum(mouse_group==1))],'.','linewidth',2,'color','b')
errorbar([4 5],[mean(mean_significant_cells_is_sig_within_env(mouse_group==2)),mean(mean_not_significant_cells_is_sig_within_env(mouse_group==2))],[std(mean_significant_cells_is_sig_within_env(mouse_group==2))./sqrt(sum(mouse_group==2)),std(mean_not_significant_cells_is_sig_within_env(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','r')
xlim([0.5 5.5])
ylim([0 0.3])
ylabel('Within env. recurrence probability')
box off
hold on
axis square
set(gca,'xtick',[1 2 4 5])
set(gca,'xticklabels',{'Dep','Ind','Dep','Ind'})
set(gca,'ytick',0:0.1:0.3)
set(gca,'fontsize',16)
axis square
box off
title('Figure S4I, inset')

% Figure S4J - Average pairwise noise correlations between cells with the
% same preferred position:
x_vec_1=1+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_2=2+0.3*(0.5-rand(1,sum(mouse_group==2)));
figure
bar([1 2],[mean(mean_same_pos_noise_correlations_per_mouse_novel(mouse_group==1)),mean(mean_same_pos_noise_correlations_per_mouse_novel(mouse_group==2))],0.5,'FaceColor','none')
hold on
errorbar([1 2],[mean(mean_same_pos_noise_correlations_per_mouse_novel(mouse_group==1)),mean(mean_same_pos_noise_correlations_per_mouse_novel(mouse_group==2))],[std(mean_same_pos_noise_correlations_per_mouse_novel(mouse_group==1))./sqrt(sum(mouse_group==1)),std(mean_same_pos_noise_correlations_per_mouse_novel(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','k')
scatter(x_vec_1,mean_same_pos_noise_correlations_per_mouse_novel(mouse_group==1),25,'k','filled');
scatter(x_vec_2,mean_same_pos_noise_correlations_per_mouse_novel(mouse_group==2),25,'k','filled');
xlim([0.5 2.5])
ylim([0 0.02])
ylabel('Noise correlation')
box off
hold on
axis square
set(gca,'xtick',[1 2])
set(gca,'xticklabels',{'CA1','CA3'})
set(gca,'fontsize',16)
axis square
box off
title('Figure S4J - Same preferred position')

% Figure S4K - Pairwise noise correlations versus difference in preferred positions:
figure
errorbar(0:bin_size:bin_size*(num_bins-1),mean_noise_corr_versus_shift_CA1_novel,std_noise_corr_versus_shift_CA1_novel./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
hold on
errorbar(0:bin_size:bin_size*(num_bins-1),mean_noise_corr_versus_shift_CA3_novel,std_noise_corr_versus_shift_CA3_novel./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
hold on
errorbar(0:bin_size:bin_size*(num_bins-1),mean_shuffle_noise_corr_versus_shift_novel,std_shuffle_noise_corr_versus_shift_novel./sqrt(sum(mouse_group==1 | mouse_group==2)),'-*k','linewidth',2)
xlim([0 40])
ylim([-0.001 0.01])
set(gca,'ytick',0:0.005:0.02)
set(gca,'yticklabels',0:0.005:0.02)
xlabel('Positional difference (cm)')
ylabel('Noise correlation')
set(gca,'fontsize',16)
legend('CA1','CA3','Shuffle')
legend('boxoff')
box off
axis square
title('Figure S4K')

% Figure S4L - Pairwise noise correlations versus distance in brain tissue:
centers_distance_vector=zeros(1,length(distance_vector)-1);
for n=1:length(distance_vector)-1
    centers_distance_vector(n)=mean(distance_vector(n:n+1));
end
figure
errorbar(centers_distance_vector(2:13),mean_noise_correlation_versus_distance_CA1_novel(2:13),std_noise_correlation_versus_distance_CA1_novel(2:13)./sqrt(sum(mouse_group==1)),'-*b','linewidth',2)
hold on
errorbar(centers_distance_vector(2:13),mean_noise_correlation_versus_distance_CA3_novel(2:13),std_noise_correlation_versus_distance_CA3_novel(2:13)./sqrt(sum(mouse_group==2)),'-*r','linewidth',2)
errorbar(centers_distance_vector(2:13),mean_shuffle_noise_correlation_versus_distance_novel(2:13),std_shuffle_noise_correlation_versus_distance_novel(2:13)./sqrt(sum(mouse_group<3)),'-*k','linewidth',2)
xlim([0 250])
ylim([-0.001 0.01])
xlabel('Distance (microns)')
ylabel('Noise correlation')
set(gca,'fontsize',16)
legend('CA1','CA3','Shuffle')
legend('boxoff')
box off
axis square
set(gca,'fontsize',16)
box off
axis square
title('Figure S4L')


%% Figure S5 - Comparison of coding properties between peer-dependent and peer-independent place cells:

% Figure S5A-C - Coding properties versus firing rates:
SI_SSR_single_cells=cell(1,number_of_mice);
average_rates_single_cells=cell(1,number_of_mice);
within_day_tuning_corr_single_cells=cell(1,number_of_mice);

for n=1:number_of_mice
    for k=1:maximal_number_of_sessions
        within_day_tuning_corr_single_cells{n}=[within_day_tuning_corr_single_cells{n},across_mice_data{n}.within_day_within_day_tuning_corr_single_cells{k}];
        SI_SSR_single_cells{n}=[SI_SSR_single_cells{n},across_mice_data{n}.within_day_SI_SSR_single_cells{k}];
        average_rates_single_cells{n}=[average_rates_single_cells{n},across_mice_data{n}.within_day_average_rates_single_cells{k}];
    end
end

for n=1:number_of_mice
    within_day_tuning_corr_single_cells{n}(isnan(average_rates_single_cells{n}))=[];
    SI_SSR_single_cells{n}(isnan(average_rates_single_cells{n}))=[];
    average_rates_single_cells{n}(isnan(average_rates_single_cells{n}))=[];
    
    within_day_tuning_corr_single_cells{n}(isnan(SI_SSR_single_cells{n}))=[];
    average_rates_single_cells{n}(isnan(SI_SSR_single_cells{n}))=[];
    SI_SSR_single_cells{n}(isnan(SI_SSR_single_cells{n}))=[];
    
    average_rates_single_cells{n}(isnan(within_day_tuning_corr_single_cells{n}))=[];
    SI_SSR_single_cells{n}(isnan(within_day_tuning_corr_single_cells{n}))=[];
    within_day_tuning_corr_single_cells{n}(isnan(within_day_tuning_corr_single_cells{n}))=[];
end

rate_vec=0:0.4:3.6;
average_rate_vec=nan(1,length(rate_vec)-1);
average_within_information_per_rate=nan(number_of_mice,length(rate_vec)-1);
average_within_correlation_per_rate=nan(number_of_mice,length(rate_vec)-1);
for k=1:number_of_mice
    for n=1:length(rate_vec)-1
        this_rate_indexes=find(average_rates_single_cells{k}>=rate_vec(n) & average_rates_single_cells{k}<=rate_vec(n+1));
        average_within_information_per_rate(k,n)=mean(SI_SSR_single_cells{k}(this_rate_indexes),'omitnan');
        average_within_correlation_per_rate(k,n)=mean(within_day_tuning_corr_single_cells{k}(this_rate_indexes),'omitnan');
        if k==1
            average_rate_vec(n)=mean(rate_vec(n:n+1));
        end
    end
end

mean_within_information_per_rate_CA1=mean(average_within_information_per_rate(mouse_group==1,:),'omitnan');
sd_within_information_per_rate_CA1=std(average_within_information_per_rate(mouse_group==1,:),'omitnan');
mean_within_information_per_rate_CA3=mean(average_within_information_per_rate(mouse_group==2,:),'omitnan');
sd_within_information_per_rate_CA3=std(average_within_information_per_rate(mouse_group==2,:),'omitnan');

mean_within_correlation_per_rate_CA1=mean(average_within_correlation_per_rate(mouse_group==1,:),'omitnan');
sd_within_correlation_per_rate_CA1=std(average_within_correlation_per_rate(mouse_group==1,:),'omitnan');
mean_within_correlation_per_rate_CA3=mean(average_within_correlation_per_rate(mouse_group==2,:),'omitnan');
sd_within_correlation_per_rate_CA3=std(average_within_correlation_per_rate(mouse_group==2,:),'omitnan');

across_days_tuning_corr_single_cells=cell(1,number_of_mice);
across_days_average_rates_single_cells=cell(1,number_of_mice);
for n=1:number_of_mice
    for k=1:maximal_number_of_sessions-1
        across_days_tuning_corr_single_cells{n}=[across_days_tuning_corr_single_cells{n},across_mice_data{n}.acorss_days_across_days_tuning_corr_single_cells{k}];
        across_days_average_rates_single_cells{n}=[across_days_average_rates_single_cells{n},across_mice_data{n}.across_days_average_rates_single_cells{k}];
    end
end

for n=1:number_of_mice
    across_days_average_rates_single_cells{n}(isnan(across_days_tuning_corr_single_cells{n}))=[];
    across_days_tuning_corr_single_cells{n}(isnan(across_days_tuning_corr_single_cells{n}))=[];
    
    across_days_tuning_corr_single_cells{n}(isnan(across_days_average_rates_single_cells{n}))=[];
    across_days_average_rates_single_cells{n}(isnan(across_days_average_rates_single_cells{n}))=[];
end

rate_vec=0:0.4:3.6;
average_across_correlation_per_rate=nan(number_of_mice,length(rate_vec)-1);
for k=1:number_of_mice
    for n=1:length(rate_vec)-1
        this_rate_indexes=find(across_days_average_rates_single_cells{k}>=rate_vec(n) & across_days_average_rates_single_cells{k}<=rate_vec(n+1));
        average_across_correlation_per_rate(k,n)=mean(across_days_tuning_corr_single_cells{k}(this_rate_indexes),'omitnan');
    end
end

mean_across_correlation_per_rate_CA1=mean(average_across_correlation_per_rate(mouse_group==1,:),'omitnan');
sd_across_correlation_per_rate_CA1=std(average_across_correlation_per_rate(mouse_group==1,:),'omitnan');
mean_across_correlation_per_rate_CA3=mean(average_across_correlation_per_rate(mouse_group==2,:),'omitnan');
sd_across_correlation_per_rate_CA3=std(average_across_correlation_per_rate(mouse_group==2,:),'omitnan');

% Figure S5D-I - Comparison of coding properties between peer-dependent
% and peer-independent cells:
average_number_of_fields_dynamics_dependent=nan(number_of_mice,maximal_number_of_sessions,num_trials);
average_max_field_size_dynamics_dependent=nan(number_of_mice,maximal_number_of_sessions,num_trials);
average_number_of_fields_dynamics_independent=nan(number_of_mice,maximal_number_of_sessions,num_trials);
average_max_field_size_dynamics_independent=nan(number_of_mice,maximal_number_of_sessions,num_trials);
for mouse_index=1:number_of_mice
    average_number_of_fields_dynamics_dependent(mouse_index,:,:)=across_mice_data{mouse_index}.average_number_of_fields_dynamics_dependent;
    average_max_field_size_dynamics_dependent(mouse_index,:,:)=across_mice_data{mouse_index}.average_max_field_size_dynamics_dependent;
    average_number_of_fields_dynamics_independent(mouse_index,:,:)=across_mice_data{mouse_index}.average_number_of_fields_dynamics_independent;
    average_max_field_size_dynamics_independent(mouse_index,:,:)=across_mice_data{mouse_index}.average_max_field_size_dynamics_independent;
end
average_number_of_fields_days_dependent=mean(average_number_of_fields_dynamics_dependent,3,'omitnan');
average_max_field_size_days_dependent=mean(average_max_field_size_dynamics_dependent,3,'omitnan');
average_number_of_fields_days_independent=mean(average_number_of_fields_dynamics_independent,3,'omitnan');
average_max_field_size_days_independent=mean(average_max_field_size_dynamics_independent,3,'omitnan');

average_number_of_fields_dependent=mean(average_number_of_fields_days_dependent,2,'omitnan');
average_max_field_size_dependent=mean(average_max_field_size_days_dependent,2,'omitnan');
average_number_of_fields_independent=mean(average_number_of_fields_days_independent,2,'omitnan');
average_max_field_size_independent=mean(average_max_field_size_days_independent,2,'omitnan');

rate_matrix_correlated_cells=nan(number_of_mice,maximal_number_of_sessions*num_trials,maximal_number_of_sessions*num_trials);
rate_matrix_uncorrelated_cells=nan(number_of_mice,maximal_number_of_sessions*num_trials,maximal_number_of_sessions*num_trials);
for n=1:number_of_mice
    rate_matrix_correlated_cells(n,:,:)=across_mice_data{n}.rate_correlations_correlated_cells;
    rate_matrix_uncorrelated_cells(n,:,:)=across_mice_data{n}.rate_correlations_uncorrelated_cells;
end

rate_consecutive_days_correlated=nan(number_of_mice,maximal_number_of_sessions);
rate_consecutive_days_uncorrelated=nan(number_of_mice,maximal_number_of_sessions);
for m=1:number_of_mice
    for n=1:maximal_number_of_sessions/2-1
        rate_consecutive_days_correlated(m,n+1)=mean([rate_matrix_correlated_cells(m,(n-1)*num_trials+1,n*num_trials+1), rate_matrix_correlated_cells(m,(n-1)*num_trials+2,n*num_trials+2)]);
        rate_consecutive_days_uncorrelated(m,n+1)=mean([rate_matrix_uncorrelated_cells(m,(n-1)*num_trials+1,n*num_trials+1), rate_matrix_uncorrelated_cells(m,(n-1)*num_trials+2,n*num_trials+2)]);
        rate_consecutive_days_correlated(m,n+maximal_number_of_sessions/2+1)=mean([rate_matrix_correlated_cells(m,(n+maximal_number_of_sessions/2-1)*num_trials+1,(n+maximal_number_of_sessions/2)*num_trials+1), rate_matrix_correlated_cells(m,(n+maximal_number_of_sessions/2-1)*num_trials+2,(n+maximal_number_of_sessions/2)*num_trials+2)]);
        rate_consecutive_days_uncorrelated(m,n+maximal_number_of_sessions/2+1)=mean([rate_matrix_uncorrelated_cells(m,(n+maximal_number_of_sessions/2-1)*num_trials+1,(n+maximal_number_of_sessions/2)*num_trials+1), rate_matrix_uncorrelated_cells(m,(n+maximal_number_of_sessions/2-1)*num_trials+2,(n+maximal_number_of_sessions/2)*num_trials+2)]);
    end
end

all_noise_correlations_same_pos_corr_cells_per_mouse=cell(1,number_of_mice);
all_noise_correlations_same_pos_uncorr_cells_per_mouse=cell(1,number_of_mice);

for n=1:number_of_mice
    for k=1:maximal_number_of_sessions
        all_noise_correlations_same_pos_corr_cells_per_mouse{n}=[all_noise_correlations_same_pos_corr_cells_per_mouse{n},across_mice_data{n}.noise_correlations_same_pos_correlated_single_cells{k}];
        all_noise_correlations_same_pos_uncorr_cells_per_mouse{n}=[all_noise_correlations_same_pos_uncorr_cells_per_mouse{n},across_mice_data{n}.noise_correlations_same_pos_uncorrelated_single_cells{k}];
    end
end

mean_noise_correlations_same_pos_corr_cells_per_mouse=nan(1,number_of_mice);
mean_noise_correlations_same_pos_other_cells_per_mouse=nan(1,number_of_mice);
for n=1:number_of_mice
    mean_noise_correlations_same_pos_corr_cells_per_mouse(n)=mean(all_noise_correlations_same_pos_corr_cells_per_mouse{n},'omitnan');
    mean_noise_correlations_same_pos_other_cells_per_mouse(n)=mean([all_noise_correlations_same_pos_uncorr_cells_per_mouse{n},all_noise_correlations_same_pos_uncorr_cells_per_mouse{n}],'omitnan');
end

average_rates_correlated_sessions=nan(number_of_mice,maximal_number_of_sessions);
average_rates_uncorrelated_sessions=nan(number_of_mice,maximal_number_of_sessions);
for n=1:number_of_mice
    average_rates_correlated_sessions(n,:)=across_mice_data{n}.average_rates_correlated_sessions;
    average_rates_uncorrelated_sessions(n,:)=across_mice_data{n}.average_rates_uncorrelated_sessions;
end

mean_per_mouse_rate_consecutive_days_correlated=mean(rate_consecutive_days_correlated,2,'omitnan');
mean_per_mouse_rate_consecutive_days_uncorrelated=mean(rate_consecutive_days_uncorrelated,2,'omitnan');
average_rates_per_mouse_correlated_cells=mean(average_rates_correlated_sessions,2,'omitnan');
average_rates_per_mouse_uncorrelated_cells=mean(average_rates_uncorrelated_sessions,2,'omitnan');

DKL_peer_dependent=nan(2,2,number_of_mice);
DKL_peer_independent=nan(2,2,number_of_mice);
for n=1:number_of_mice
    this_mouse_dependent_dist=hist(across_mice_data{n}.all_correlated_preferred_positions_right{1},1:num_bins);
    this_mouse_dependent_dist=this_mouse_dependent_dist./sum(this_mouse_dependent_dist);
    DKL_peer_dependent(1,1,n)=-sum(this_mouse_dependent_dist.*log2(1/num_bins*ones(1,num_bins)./this_mouse_dependent_dist));
    this_mouse_dependent_dist=hist(across_mice_data{n}.all_correlated_preferred_positions_left{1},1:num_bins);
    this_mouse_dependent_dist=this_mouse_dependent_dist./sum(this_mouse_dependent_dist);
    DKL_peer_dependent(1,2,n)=-sum(this_mouse_dependent_dist.*log2(1/num_bins*ones(1,num_bins)./this_mouse_dependent_dist));
    
    this_mouse_independent_dist=hist(across_mice_data{n}.all_uncorrelated_preferred_positions_right{1},1:num_bins);
    this_mouse_independent_dist=this_mouse_independent_dist./sum(this_mouse_independent_dist);
    DKL_peer_independent(1,1,n)=-sum(this_mouse_independent_dist.*log2(1/num_bins*ones(1,num_bins)./this_mouse_independent_dist));
    this_mouse_independent_dist=hist(across_mice_data{n}.all_uncorrelated_preferred_positions_left{1},1:num_bins);
    this_mouse_independent_dist=this_mouse_independent_dist./sum(this_mouse_independent_dist);
    DKL_peer_independent(1,2,n)=-sum(this_mouse_independent_dist.*log2(1/num_bins*ones(1,num_bins)./this_mouse_independent_dist));
    
    
    this_mouse_dependent_dist=hist(across_mice_data{n}.all_correlated_preferred_positions_right{2},1:num_bins);
    this_mouse_dependent_dist=this_mouse_dependent_dist./sum(this_mouse_dependent_dist);
    DKL_peer_dependent(2,1,n)=-sum(this_mouse_dependent_dist.*log2(1/num_bins*ones(1,num_bins)./this_mouse_dependent_dist));
    this_mouse_dependent_dist=hist(across_mice_data{n}.all_correlated_preferred_positions_left{2},1:num_bins);
    this_mouse_dependent_dist=this_mouse_dependent_dist./sum(this_mouse_dependent_dist);
    DKL_peer_dependent(2,2,n)=-sum(this_mouse_dependent_dist.*log2(1/num_bins*ones(1,num_bins)./this_mouse_dependent_dist));
    
    this_mouse_independent_dist=hist(across_mice_data{n}.all_uncorrelated_preferred_positions_right{2},1:num_bins);
    this_mouse_independent_dist=this_mouse_independent_dist./sum(this_mouse_independent_dist);
    DKL_peer_independent(2,1,n)=-sum(this_mouse_independent_dist.*log2(1/num_bins*ones(1,num_bins)./this_mouse_independent_dist));
    this_mouse_independent_dist=hist(across_mice_data{n}.all_uncorrelated_preferred_positions_left{2},1:num_bins);
    this_mouse_independent_dist=this_mouse_independent_dist./sum(this_mouse_independent_dist);
    DKL_peer_independent(2,2,n)=-sum(this_mouse_independent_dist.*log2(1/num_bins*ones(1,num_bins)./this_mouse_independent_dist));
end
average_DKL_peer_dependent=mean(mean(DKL_peer_dependent));
average_DKL_peer_independent=mean(mean(DKL_peer_independent));

% Figure S5A - Spatial information versus firing rates:
figure
errorbar(average_rate_vec,mean_within_information_per_rate_CA1,sd_within_information_per_rate_CA1./sqrt(sum(mouse_group==1)),'-*','markersize',5,'color','b','linewidth',2)
hold on
errorbar(average_rate_vec,mean_within_information_per_rate_CA3,sd_within_information_per_rate_CA3./sqrt(sum(mouse_group==2)),'-*','markersize',5,'color','r','linewidth',2)
xlim([0 3.6])
ylim([0 2.5])
xlabel('Firing rate (spike/sec)')
ylabel('Spatial information (bit/spike)')
legend('CA1','CA3','Location','Northwest')
legend boxoff
box off
axis square
set(gca,'fontsize',16)
title('Figure S5A')

% Figure S5B - Within day tuning curve correlation versus firing rates:
figure
errorbar(average_rate_vec,mean_within_correlation_per_rate_CA1,sd_within_correlation_per_rate_CA1./sqrt(sum(mouse_group==1)),'-*','markersize',5,'color','b','linewidth',2)
hold on
errorbar(average_rate_vec,mean_within_correlation_per_rate_CA3,sd_within_correlation_per_rate_CA3./sqrt(sum(mouse_group==2)),'-*','markersize',5,'color','r','linewidth',2)
xlim([0 3.6])
ylim([0 1])
xlabel('Firing rate (spike/sec)')
ylabel('Tuning curve correlation')
legend('CA1','CA3','Location','Southwest')
legend boxoff
box off
axis square
set(gca,'fontsize',16)
title('Figure S5B - Within day')

% Figure S5C - Across days tuning curve correlation versus firing rates:
figure
errorbar(average_rate_vec,mean_across_correlation_per_rate_CA1,sd_across_correlation_per_rate_CA1./sqrt(sum(mouse_group==1)),'-*','markersize',5,'color','b','linewidth',2)
hold on
errorbar(average_rate_vec,mean_across_correlation_per_rate_CA3,sd_across_correlation_per_rate_CA3./sqrt(sum(mouse_group==2)),'-*','markersize',5,'color','r','linewidth',2)
xlim([0 3.6])
ylim([0 1])
xlabel('Firing rate (spike/sec)')
ylabel('Tuning curve correlation')
legend('CA1','CA3','Location','Northwest')
legend boxoff
box off
axis square
set(gca,'fontsize',16)
title('Figure S5C - Across days')

% Figure S5D - Average place field size for peer-dependent and peer-independent cells:
figure
x_vec_1=1+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_2=2+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_3=4+0.3*(0.5-rand(1,sum(mouse_group==2)));
x_vec_4=5+0.3*(0.5-rand(1,sum(mouse_group==2)));
plot([1 2],bin_size*[mean(average_max_field_size_dependent(mouse_group==1)), mean(average_max_field_size_independent(mouse_group==1))],'.-b','markersize',30,'linewidth',2)
hold on
plot([4 5],bin_size*[mean(average_max_field_size_dependent(mouse_group==2)), mean(average_max_field_size_independent(mouse_group==2))],'.-r','markersize',30,'linewidth',2)
errorbar([1 2],bin_size*[mean(average_max_field_size_dependent(mouse_group==1)),mean(average_max_field_size_independent(mouse_group==1))],bin_size*[std(average_max_field_size_dependent(mouse_group==1))./sqrt(sum(mouse_group==1)),std(average_max_field_size_independent(mouse_group==1))./sqrt(sum(mouse_group==1))],'.','linewidth',2,'color','b')
errorbar([4 5],bin_size*[mean(average_max_field_size_dependent(mouse_group==2)),mean(average_max_field_size_independent(mouse_group==2))],bin_size*[std(average_max_field_size_dependent(mouse_group==2))./sqrt(sum(mouse_group==2)),std(average_max_field_size_independent(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','r')
for n=1:sum(mouse_group==1)
    plot([x_vec_1(n), x_vec_2(n)],bin_size*[average_max_field_size_dependent(CA1_indexes(n)),average_max_field_size_independent(CA1_indexes(n))],'-','color',[0.7 0.7 0.7])
end
for n=1:sum(mouse_group==2)
    plot([x_vec_3(n), x_vec_4(n)],bin_size*[average_max_field_size_dependent(CA3_indexes(n)),average_max_field_size_independent(CA3_indexes(n))],'-','color',[0.7 0.7 0.7])
end
plot([1 2],bin_size*[mean(average_max_field_size_dependent(mouse_group==1)), mean(average_max_field_size_independent(mouse_group==1))],'.-b','markersize',30,'linewidth',2)
plot([4 5],bin_size*[mean(average_max_field_size_dependent(mouse_group==2)), mean(average_max_field_size_independent(mouse_group==2))],'.-r','markersize',30,'linewidth',2)
errorbar([1 2],bin_size*[mean(average_max_field_size_dependent(mouse_group==1)),mean(average_max_field_size_independent(mouse_group==1))],bin_size*[std(average_max_field_size_dependent(mouse_group==1))./sqrt(sum(mouse_group==1)),std(average_max_field_size_independent(mouse_group==1))./sqrt(sum(mouse_group==1))],'.','linewidth',2,'color','b')
errorbar([4 5],bin_size*[mean(average_max_field_size_dependent(mouse_group==2)),mean(average_max_field_size_independent(mouse_group==2))],bin_size*[std(average_max_field_size_dependent(mouse_group==2))./sqrt(sum(mouse_group==2)),std(average_max_field_size_independent(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','r')
xlim([0.5 5.5])
ylim([0 15])
ylabel('Place field size (cm)')
box off
hold on
axis square
set(gca,'xtick',[1 2 4 5])
set(gca,'xticklabels',{'Dep','Ind','Dep','Ind'})
set(gca,'fontsize',16)
axis square
box off
title('Figure S5D')

% Figure S5E - Average number of place fields per cell peer-dependent and peer-independent cells:
figure
x_vec_1=1+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_2=2+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_3=4+0.3*(0.5-rand(1,sum(mouse_group==2)));
x_vec_4=5+0.3*(0.5-rand(1,sum(mouse_group==2)));
plot([1 2],[mean(average_number_of_fields_dependent(mouse_group==1)), mean(average_number_of_fields_independent(mouse_group==1))],'.-b','markersize',30,'linewidth',2)
hold on
plot([4 5],[mean(average_number_of_fields_dependent(mouse_group==2)), mean(average_number_of_fields_independent(mouse_group==2))],'.-r','markersize',30,'linewidth',2)
errorbar([1 2],[mean(average_number_of_fields_dependent(mouse_group==1)),mean(average_number_of_fields_independent(mouse_group==1))],[std(average_number_of_fields_dependent(mouse_group==1))./sqrt(sum(mouse_group==1)),std(average_number_of_fields_independent(mouse_group==1))./sqrt(sum(mouse_group==1))],'.','linewidth',2,'color','b')
errorbar([4 5],[mean(average_number_of_fields_dependent(mouse_group==2)),mean(average_number_of_fields_independent(mouse_group==2))],[std(average_number_of_fields_dependent(mouse_group==2))./sqrt(sum(mouse_group==2)),std(average_number_of_fields_independent(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','r')
for n=1:sum(mouse_group==1)
    plot([x_vec_1(n), x_vec_2(n)],[average_number_of_fields_dependent(CA1_indexes(n)),average_number_of_fields_independent(CA1_indexes(n))],'-','color',[0.7 0.7 0.7])
end
for n=1:sum(mouse_group==2)
    plot([x_vec_3(n), x_vec_4(n)],[average_number_of_fields_dependent(CA3_indexes(n)),average_number_of_fields_independent(CA3_indexes(n))],'-','color',[0.7 0.7 0.7])
end
plot([1 2],[mean(average_number_of_fields_dependent(mouse_group==1)), mean(average_number_of_fields_independent(mouse_group==1))],'.-b','markersize',30,'linewidth',2)
plot([4 5],[mean(average_number_of_fields_dependent(mouse_group==2)), mean(average_number_of_fields_independent(mouse_group==2))],'.-r','markersize',30,'linewidth',2)
errorbar([1 2],[mean(average_number_of_fields_dependent(mouse_group==1)),mean(average_number_of_fields_independent(mouse_group==1))],[std(average_number_of_fields_dependent(mouse_group==1))./sqrt(sum(mouse_group==1)),std(average_number_of_fields_independent(mouse_group==1))./sqrt(sum(mouse_group==1))],'.','linewidth',2,'color','b')
errorbar([4 5],[mean(average_number_of_fields_dependent(mouse_group==2)),mean(average_number_of_fields_independent(mouse_group==2))],[std(average_number_of_fields_dependent(mouse_group==2))./sqrt(sum(mouse_group==2)),std(average_number_of_fields_independent(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','r')
xlim([0.5 5.5])
ylim([0 3])
ylabel('Number of fields per cell')
box off
hold on
axis square
set(gca,'xtick',[1 2 4 5])
set(gca,'xticklabels',{'Dep','Ind','Dep','Ind'})
set(gca,'ytick',0:1:3)
set(gca,'fontsize',16)
axis square
box off
title('Figure S5E')

% Figure S5F - Average rate correlations across days for peer-dependent and peer-independent cells:
figure
x_vec_1=1+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_2=2+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_3=4+0.3*(0.5-rand(1,sum(mouse_group==2)));
x_vec_4=5+0.3*(0.5-rand(1,sum(mouse_group==2)));
plot([1 2],[mean(mean_per_mouse_rate_consecutive_days_correlated(mouse_group==1)), mean(mean_per_mouse_rate_consecutive_days_uncorrelated(mouse_group==1))],'.-b','markersize',30,'linewidth',2)
hold on
plot([4 5],[mean(mean_per_mouse_rate_consecutive_days_correlated(mouse_group==2)), mean(mean_per_mouse_rate_consecutive_days_uncorrelated(mouse_group==2))],'.-r','markersize',30,'linewidth',2)
errorbar([1 2],[mean(mean_per_mouse_rate_consecutive_days_correlated(mouse_group==1)),mean(mean_per_mouse_rate_consecutive_days_uncorrelated(mouse_group==1))],[std(mean_per_mouse_rate_consecutive_days_correlated(mouse_group==1))./sqrt(sum(mouse_group==1)),std(mean_per_mouse_rate_consecutive_days_uncorrelated(mouse_group==1))./sqrt(sum(mouse_group==1))],'.','linewidth',2,'color','b')
errorbar([4 5],[mean(mean_per_mouse_rate_consecutive_days_correlated(mouse_group==2)),mean(mean_per_mouse_rate_consecutive_days_uncorrelated(mouse_group==2))],[std(mean_per_mouse_rate_consecutive_days_correlated(mouse_group==2))./sqrt(sum(mouse_group==2)),std(mean_per_mouse_rate_consecutive_days_uncorrelated(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','r')
for n=1:sum(mouse_group==1)
    plot([x_vec_1(n), x_vec_2(n)],[mean_per_mouse_rate_consecutive_days_correlated(CA1_indexes(n)),mean_per_mouse_rate_consecutive_days_uncorrelated(CA1_indexes(n))],'-','color',[0.7 0.7 0.7])
end
for n=1:sum(mouse_group==2)
    plot([x_vec_3(n), x_vec_4(n)],[mean_per_mouse_rate_consecutive_days_correlated(CA3_indexes(n)),mean_per_mouse_rate_consecutive_days_uncorrelated(CA3_indexes(n))],'-','color',[0.7 0.7 0.7])
end
plot([1 2],[mean(mean_per_mouse_rate_consecutive_days_correlated(mouse_group==1)), mean(mean_per_mouse_rate_consecutive_days_uncorrelated(mouse_group==1))],'.-b','markersize',30,'linewidth',2)
plot([4 5],[mean(mean_per_mouse_rate_consecutive_days_correlated(mouse_group==2)), mean(mean_per_mouse_rate_consecutive_days_uncorrelated(mouse_group==2))],'.-r','markersize',30,'linewidth',2)
errorbar([1 2],[mean(mean_per_mouse_rate_consecutive_days_correlated(mouse_group==1)),mean(mean_per_mouse_rate_consecutive_days_uncorrelated(mouse_group==1))],[std(mean_per_mouse_rate_consecutive_days_correlated(mouse_group==1))./sqrt(sum(mouse_group==1)),std(mean_per_mouse_rate_consecutive_days_uncorrelated(mouse_group==1))./sqrt(sum(mouse_group==1))],'.','linewidth',2,'color','b')
errorbar([4 5],[mean(mean_per_mouse_rate_consecutive_days_correlated(mouse_group==2)),mean(mean_per_mouse_rate_consecutive_days_uncorrelated(mouse_group==2))],[std(mean_per_mouse_rate_consecutive_days_correlated(mouse_group==2))./sqrt(sum(mouse_group==2)),std(mean_per_mouse_rate_consecutive_days_uncorrelated(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','r')
xlim([0.5 5.5])
ylim([0 1])
ylabel('Ensemble ate correlation')
box off
hold on
axis square
set(gca,'xtick',[1 2 4 5])
set(gca,'xticklabels',{'Dep','Ind','Dep','Ind'})
set(gca,'fontsize',16)
axis square
box off
title('Figure S5F')

% Figure S5G - Average pairwise noise correlations for peer-dependent and peer-independent cells:
figure
x_vec_1=1+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_2=2+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_3=4+0.3*(0.5-rand(1,sum(mouse_group==2)));
x_vec_4=5+0.3*(0.5-rand(1,sum(mouse_group==2)));
plot([1 2],[mean(mean_noise_correlations_same_pos_corr_cells_per_mouse(mouse_group==1)), mean(mean_noise_correlations_same_pos_other_cells_per_mouse(mouse_group==1))],'.-b','markersize',30,'linewidth',2)
hold on
plot([4 5],[mean(mean_noise_correlations_same_pos_corr_cells_per_mouse(mouse_group==2)), mean(mean_noise_correlations_same_pos_other_cells_per_mouse(mouse_group==2))],'.-r','markersize',30,'linewidth',2)
errorbar([1 2],[mean(mean_noise_correlations_same_pos_corr_cells_per_mouse(mouse_group==1)),mean(mean_noise_correlations_same_pos_other_cells_per_mouse(mouse_group==1))],[std(mean_noise_correlations_same_pos_corr_cells_per_mouse(mouse_group==1))./sqrt(sum(mouse_group==1)),std(mean_noise_correlations_same_pos_other_cells_per_mouse(mouse_group==1))./sqrt(sum(mouse_group==1))],'.','linewidth',2,'color','b')
errorbar([4 5],[mean(mean_noise_correlations_same_pos_corr_cells_per_mouse(mouse_group==2)),mean(mean_noise_correlations_same_pos_other_cells_per_mouse(mouse_group==2))],[std(mean_noise_correlations_same_pos_corr_cells_per_mouse(mouse_group==2))./sqrt(sum(mouse_group==2)),std(mean_noise_correlations_same_pos_other_cells_per_mouse(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','r')
for n=1:sum(mouse_group==1)
    plot([x_vec_1(n), x_vec_2(n)],[mean_noise_correlations_same_pos_corr_cells_per_mouse(CA1_indexes(n)),mean_noise_correlations_same_pos_other_cells_per_mouse(CA1_indexes(n))],'-','color',[0.7 0.7 0.7])
end
for n=1:sum(mouse_group==2)
    plot([x_vec_3(n), x_vec_4(n)],[mean_noise_correlations_same_pos_corr_cells_per_mouse(CA3_indexes(n)),mean_noise_correlations_same_pos_other_cells_per_mouse(CA3_indexes(n))],'-','color',[0.7 0.7 0.7])
end
plot([1 2],[mean(mean_noise_correlations_same_pos_corr_cells_per_mouse(mouse_group==1)), mean(mean_noise_correlations_same_pos_other_cells_per_mouse(mouse_group==1))],'.-b','markersize',30,'linewidth',2)
plot([4 5],[mean(mean_noise_correlations_same_pos_corr_cells_per_mouse(mouse_group==2)), mean(mean_noise_correlations_same_pos_other_cells_per_mouse(mouse_group==2))],'.-r','markersize',30,'linewidth',2)
errorbar([1 2],[mean(mean_noise_correlations_same_pos_corr_cells_per_mouse(mouse_group==1)),mean(mean_noise_correlations_same_pos_other_cells_per_mouse(mouse_group==1))],[std(mean_noise_correlations_same_pos_corr_cells_per_mouse(mouse_group==1))./sqrt(sum(mouse_group==1)),std(mean_noise_correlations_same_pos_other_cells_per_mouse(mouse_group==1))./sqrt(sum(mouse_group==1))],'.','linewidth',2,'color','b')
errorbar([4 5],[mean(mean_noise_correlations_same_pos_corr_cells_per_mouse(mouse_group==2)),mean(mean_noise_correlations_same_pos_other_cells_per_mouse(mouse_group==2))],[std(mean_noise_correlations_same_pos_corr_cells_per_mouse(mouse_group==2))./sqrt(sum(mouse_group==2)),std(mean_noise_correlations_same_pos_other_cells_per_mouse(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','r')
xlim([0.5 5.5])
ylim([0 0.03])
ylabel('Noise correlation')
box off
hold on
axis square
set(gca,'ytick',0:0.01:0.03)
set(gca,'xtick',[1 2 4 5])
set(gca,'xticklabels',{'Dep','Ind','Dep','Ind'})
set(gca,'fontsize',16)
axis square
box off
title('Figure S5G')

% Figure S5H - Average firing rates for peer-dependent and peer-independent cells:
figure
x_vec_1=1+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_2=2+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_3=4+0.3*(0.5-rand(1,sum(mouse_group==2)));
x_vec_4=5+0.3*(0.5-rand(1,sum(mouse_group==2)));
plot([1 2],[mean(average_rates_per_mouse_correlated_cells(mouse_group==1)), mean(average_rates_per_mouse_uncorrelated_cells(mouse_group==1))],'.-b','markersize',30,'linewidth',2)
hold on
plot([4 5],[mean(average_rates_per_mouse_correlated_cells(mouse_group==2)), mean(average_rates_per_mouse_uncorrelated_cells(mouse_group==2))],'.-r','markersize',30,'linewidth',2)
errorbar([1 2],[mean(average_rates_per_mouse_correlated_cells(mouse_group==1)),mean(average_rates_per_mouse_uncorrelated_cells(mouse_group==1))],[std(average_rates_per_mouse_correlated_cells(mouse_group==1))./sqrt(sum(mouse_group==1)),std(average_rates_per_mouse_uncorrelated_cells(mouse_group==1))./sqrt(sum(mouse_group==1))],'.','linewidth',2,'color','b')
errorbar([4 5],[mean(average_rates_per_mouse_correlated_cells(mouse_group==2)),mean(average_rates_per_mouse_uncorrelated_cells(mouse_group==2))],[std(average_rates_per_mouse_correlated_cells(mouse_group==2))./sqrt(sum(mouse_group==2)),std(average_rates_per_mouse_uncorrelated_cells(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','r')
for n=1:sum(mouse_group==1)
    plot([x_vec_1(n), x_vec_2(n)],[average_rates_per_mouse_correlated_cells(CA1_indexes(n)),average_rates_per_mouse_uncorrelated_cells(CA1_indexes(n))],'-','color',[0.7 0.7 0.7])
end
for n=1:sum(mouse_group==2)
    plot([x_vec_3(n), x_vec_4(n)],[average_rates_per_mouse_correlated_cells(CA3_indexes(n)),average_rates_per_mouse_uncorrelated_cells(CA3_indexes(n))],'-','color',[0.7 0.7 0.7])
end
plot([1 2],[mean(average_rates_per_mouse_correlated_cells(mouse_group==1)), mean(average_rates_per_mouse_uncorrelated_cells(mouse_group==1))],'.-b','markersize',30,'linewidth',2)
plot([4 5],[mean(average_rates_per_mouse_correlated_cells(mouse_group==2)), mean(average_rates_per_mouse_uncorrelated_cells(mouse_group==2))],'.-r','markersize',30,'linewidth',2)
errorbar([1 2],[mean(average_rates_per_mouse_correlated_cells(mouse_group==1)),mean(average_rates_per_mouse_uncorrelated_cells(mouse_group==1))],[std(average_rates_per_mouse_correlated_cells(mouse_group==1))./sqrt(sum(mouse_group==1)),std(average_rates_per_mouse_uncorrelated_cells(mouse_group==1))./sqrt(sum(mouse_group==1))],'.','linewidth',2,'color','b')
errorbar([4 5],[mean(average_rates_per_mouse_correlated_cells(mouse_group==2)),mean(average_rates_per_mouse_uncorrelated_cells(mouse_group==2))],[std(average_rates_per_mouse_correlated_cells(mouse_group==2))./sqrt(sum(mouse_group==2)),std(average_rates_per_mouse_uncorrelated_cells(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','r')
xlim([0.5 5.5])
ylim([0 2])
ylabel('Average rate (spike/sec)')
box off
hold on
axis square
set(gca,'xtick',[1 2 4 5])
set(gca,'xticklabels',{'Dep','Ind','Dep','Ind'})
set(gca,'fontsize',16)
axis square
box off
title('Figure S5H')

% Figure S5I -Divergence from a uniform distribution for peer-dependent and peer-independent cells:
figure
x_vec_1=1+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_2=2+0.3*(0.5-rand(1,sum(mouse_group==1)));
x_vec_3=4+0.3*(0.5-rand(1,sum(mouse_group==2)));
x_vec_4=5+0.3*(0.5-rand(1,sum(mouse_group==2)));
plot([1 2],[mean(average_DKL_peer_dependent(mouse_group==1)), mean(average_DKL_peer_independent(mouse_group==1))],'.-b','markersize',30,'linewidth',2)
hold on
plot([4 5],[mean(average_DKL_peer_dependent(mouse_group==2)), mean(average_DKL_peer_independent(mouse_group==2))],'.-r','markersize',30,'linewidth',2)
errorbar([1 2],[mean(average_DKL_peer_dependent(mouse_group==1)),mean(average_DKL_peer_independent(mouse_group==1))],[std(average_DKL_peer_dependent(mouse_group==1))./sqrt(sum(mouse_group==1)),std(average_DKL_peer_independent(mouse_group==1))./sqrt(sum(mouse_group==1))],'.','linewidth',2,'color','b')
errorbar([4 5],[mean(average_DKL_peer_dependent(mouse_group==2)),mean(average_DKL_peer_independent(mouse_group==2))],[std(average_DKL_peer_dependent(mouse_group==2))./sqrt(sum(mouse_group==2)),std(average_DKL_peer_independent(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','r')
for n=1:sum(mouse_group==1)
    plot([x_vec_1(n), x_vec_2(n)],[average_DKL_peer_dependent(CA1_indexes(n)),average_DKL_peer_independent(CA1_indexes(n))],'-','color',[0.7 0.7 0.7])
end
for n=1:sum(mouse_group==2)
    plot([x_vec_3(n), x_vec_4(n)],[average_DKL_peer_dependent(CA3_indexes(n)),average_DKL_peer_independent(CA3_indexes(n))],'-','color',[0.7 0.7 0.7])
end
plot([1 2],[mean(average_DKL_peer_dependent(mouse_group==1)), mean(average_DKL_peer_independent(mouse_group==1))],'.-b','markersize',30,'linewidth',2)
plot([4 5],[mean(average_DKL_peer_dependent(mouse_group==2)), mean(average_DKL_peer_independent(mouse_group==2))],'.-r','markersize',30,'linewidth',2)
errorbar([1 2],[mean(average_DKL_peer_dependent(mouse_group==1)),mean(average_DKL_peer_independent(mouse_group==1))],[std(average_DKL_peer_dependent(mouse_group==1))./sqrt(sum(mouse_group==1)),std(average_DKL_peer_independent(mouse_group==1))./sqrt(sum(mouse_group==1))],'.','linewidth',2,'color','b')
errorbar([4 5],[mean(average_DKL_peer_dependent(mouse_group==2)),mean(average_DKL_peer_independent(mouse_group==2))],[std(average_DKL_peer_dependent(mouse_group==2))./sqrt(sum(mouse_group==2)),std(average_DKL_peer_independent(mouse_group==2))./sqrt(sum(mouse_group==2))],'.','linewidth',2,'color','r')
xlim([0.5 5.5])
ylim([0 0.3])
ylabel('Divergence from uniform')
box off
hold on
axis square
set(gca,'xtick',[1 2 4 5])
set(gca,'fontsize',16)
axis square
box off
title('Figure S5I')
