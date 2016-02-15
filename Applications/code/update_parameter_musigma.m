% count_completed_core = 0;
% for i = 1 : num_seq 
%     fN = [resultsFolder coreNAME '_input' num2str(i) '_bb' num2str(bb) '_updateD'];
%     if exist(fN, 'file')
%         count_completed_core = count_completed_core + 1;
%         completed_cores(count_completed_core) = i;
%     else
%         disp(['No result from ' num2str(i)])
%     end
% end

% tolSIGMA = 1e-4;
tolSIGMA = 0.025;

input_count = zeros(size(total_target, 1), num_seq);

aligned_record = cell(size(total_target, 1),1);
for musigmaloc = 1 : size(total_target, 1)
    aligned_sum = 0;
    aligned_num = 0;
    aligned_squared_sum = 0;
    
    aligned_tmp1 = [];
    aligned_tmp2 = [];
    aligned_tmp3 = [];
    
    for i = completed_cores
        % load bin_mean_index, sample_seq, and input_scaled_new
        load([resultsFolder coreNAME '_input' num2str(i) '_iter' num2str(bb) '_updateD'])
        input_scaled_new(:,2) = input_scaled_new(:,2) + meanShift;
        
        binI = find(bin_mean_index == musigmaloc);
        
        for muI = 1 : length(binI)
            [seqI, seqJ] = find(sample_seq == binI(muI));
            if isempty(seqI)
            else
                aligned_sum = aligned_sum + sum(input_scaled_new(seqI,2));
                aligned_num = aligned_num + length(seqI);
                
                aligned_tmp1 = [aligned_tmp1; i*ones(length(seqI),1)];
                aligned_tmp2 = [aligned_tmp2; seqI];
                aligned_tmp3 = [aligned_tmp3; input_scaled_new(seqI,2)];
            end
        end
        input_count(musigmaloc, i) = aligned_num;
    end
    
    aligned_record{musigmaloc} = [aligned_tmp1 aligned_tmp2 aligned_tmp3];
    
    
    if aligned_num ~= 0
        total_mu(musigmaloc) = aligned_sum/aligned_num;
    end
    
    
    for i = completed_cores
        % load bin_mean_index, sample_seq, and input_scaled_new
        load([resultsFolder coreNAME '_input' num2str(i) '_iter' num2str(bb) '_updateD'])
        input_scaled_new(:,2) = input_scaled_new(:,2) + meanShift;
        
        binI = find(bin_mean_index == musigmaloc);
        
        for muI = 1 : length(binI)
            [seqI, seqJ] = find(sample_seq == binI(muI));
            if isempty(seqI)
            else
                tempSquared = input_scaled_new(seqI,2)-total_mu(musigmaloc);
                aligned_squared_sum = aligned_squared_sum + tempSquared'*tempSquared;
            end
        end
    end
    
    if aligned_num ~= 0
        updateMS = sqrt(aligned_squared_sum/aligned_num);
        if updateMS < tolSIGMA
            fprintf('total_sigma(%d) is less than %f in iteration %d\n', musigmaloc, tolSIGMA, bb);
            total_sigma(musigmaloc) = tolSIGMA;
        else
            total_sigma(musigmaloc) = updateMS;
        end
    end
end
input_count = [input_count(:,1) diff(input_count,1,2)];

if strcmp(target_est, 'linear')
    input_count_ratio = zeros(size(total_target, 1), num_seq);
    for musigmaloc = 1 : size(total_target, 1)
        aligned_sum = 0;
        aligned_num = 0;
        
        for i = completed_cores
            % load bin_mean_index, sample_seq, and input_scaled_new
            load([resultsFolder coreNAME '_input' num2str(i) '_iter' num2str(bb) '_updateD'])
            
            binI_p = find(bin_mean_index_p(:,1) == musigmaloc);
            
            for muI = 1 : length(binI_p)
                [seqI, ~] = find(sample_seq == binI_p(muI));
                if isempty(seqI)
                else
                    tmp1 = bin_mean_index_p(binI_p(muI),2);
                    tmp2 = bin_mean_index_n(binI_p(muI),2);
                    aligned_sum = aligned_sum + sum(input_scaled_new(seqI,2))*tmp1 - length(seqI)*tmp1*tmp2*total_mu(musigmaloc+1);
                    aligned_num = aligned_num + length(seqI)*tmp1*tmp1;
                end
            end
            
            binI_n = find(bin_mean_index_n(:,1) == musigmaloc);
            
            for muI = 1 : length(binI_n)
                [seqI, ~] = find(sample_seq == binI_n(muI));
                if isempty(seqI)
                else
                    tmp1 = bin_mean_index_p(binI_n(muI),2);
                    tmp2 = bin_mean_index_n(binI_n(muI),2);
                    aligned_sum = aligned_sum + sum(input_scaled_new(seqI,2))*tmp2 - length(seqI)*tmp1*tmp2*total_mu(musigmaloc-1);
                    aligned_num = aligned_num + length(seqI)*tmp2*tmp2;
                end
            end
            input_count_ratio(musigmaloc, i) = aligned_num;
        end
        
        
        if aligned_num ~= 0
            total_mu(musigmaloc) = aligned_sum/aligned_num;
        end
    end
    input_count_ratio = [input_count_ratio(:,1) diff(input_count_ratio,1,2)];
end








































