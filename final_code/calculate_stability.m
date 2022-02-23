function [resilience, final_diff, absolute_pre, prop_plasmid] = calculate_stability(D, A_pre, A_dur, A_post, t_post, t_dur, stability_type)
resilience = zeros(D,1);
final_diff = zeros(D,1);

absolute_pre = zeros(D,1);
prop_plasmid = zeros(D,1);
store_pre = zeros(size(A_pre,1),D);
store_post = zeros(size(A_post,1),D);
store_dur = zeros(size(A_dur,1),D);

for i=1:D
    total_pre = A_pre(:,i)+A_pre(:,i+D)+A_pre(:,i+2*D)+A_pre(:,i+3*D);
    store_pre(:,i) = total_pre;
    
    total_dur = A_dur(:,i)+A_dur(:,i+D)+A_dur(:,i+2*D)+A_dur(:,i+3*D);
    store_dur(:,i) = total_dur;
    
    total_post = A_post(:,i)+A_post(:,i+D)+A_post(:,i+2*D)+A_post(:,i+3*D);
    store_post(:,i) = total_post;
    
    plasmid_freq = (A_pre(:,i+D)+A_pre(:,i+3*D))./total_pre;
    prop_plasmid(i) = plasmid_freq(end);
    
    absolute_pre(i) = total_pre(end);
    final_pre = total_pre(end);
    
    if stability_type == 1
        tmp = log10(total_dur(end)/final_pre);
        tmp(tmp>0)=0; % to stop increases counting
        [a, b] = min(tmp);
        resilience(i) = a;
    elseif stability_type == 2
        try
            tmp = [[total_post;0],[0;total_post]];
            tmp = tmp(:,1)-tmp(:,2);
            a = find(abs(tmp)<0.05);
            resilience(i) = t_post(a(1)) - t_dur(end);
        catch
            resilience(i) = t_post(end) - t_dur(end);
        end
    elseif stability_type == 3
        try
            tmp = (total_post - final_pre)/final_pre;
            a = find(abs(tmp)<0.05);
            resilience(i) = t_post(a(1)) - t_dur(end);
        catch
            resilience(i) = t_post(end) - t_dur(end);
        end
    end
    
    tmp = total_post - final_pre;
    final_diff(i) = tmp(end);
end

if stability_type == 4
    
    save_pre = store_pre(end,:);%./sum(store_pre(end,:));
    save_dur = store_dur(end,:);%./sum( store_dur(end,:));
    
   %[simps_pre_all, invsimps_pre_all] = simpson_di(save_pre);
   %[simps_dur_all, invsimps_dur_all] = simpson_di(save_dur);
   %resilience(1) = invsimps_pre_all - invsimps_dur_all;
    
    mydis = f_dis([save_pre; save_dur], 'bc');
    resilience(1) = mydis(1,2);

%   [simps_pre_subset, invsimps_pre_subset] = simpson_di(save_pre(2:end));
%   [simps_dur_subset, invsimps_dur_subset] = simpson_di(save_dur(2:end));
   %resilience(2) = invsimps_pre_subset - invsimps_dur_subset;
   
    mydis = f_dis([save_pre(2:end); save_dur(2:end)], 'bc');
    resilience(2) = mydis(1,2);
   
end



end