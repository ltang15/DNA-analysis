%{DNA analysis}%
%Description: Analyzing a segment of DNA provided by a file of numerical values of bases.
%Calculating the total protein segments, number of bases (average, min and max length).
%Calculating percentage of DNA is used for protein process. 
%Determining which is the most frequently and the least frequently used stop codons.  

clc, clear all, close all;

load('chr1_sect.mat');

n = length (dna); % total number of bases in the given DNA file
start = 0;% initialize the start point

%initialize counting variables for corresponding stop codons
tag = 0;
tga = 0; 
taa = 0;

arr = [];%initialize the array containing the length of each protein segment
stop_freq =[];%initialize the array containing the couting values of each stop codon

%% Looping all the bases in the DNA array
for k =1:3:n-2
    % when the start point is not set, so we look for the starting codon and then assign value for it 
    if (start ==0)
        if(dna(k) == 1 && dna(k+1) ==4 && dna(k+2) ==3) %starting codon ATG (143)
            start = k;
        end
    % when the start point is set, so we look for the stopping codon and then calculate the length    
    else
        %stop codons: TAG (413) or TGA (431) or TAA (411)
        if((dna(k) == 4 && dna(k+1) == 1 && dna(k+2) == 3)||(dna(k) == 4 && dna(k+1) == 3 && dna(k+2) == 1)||(dna(k) == 4 && dna(k+1) == 1 && dna(k+2) == 1))
           
            len = k - start + 3; % calculating the number of base of each protein segment
            arr = [arr, len]; %add the length element to the array
            start = 0; % reset the start point after calculating
            
            % couting three types of stop codons
            if(dna(k) == 4 && dna(k+1) == 1 && dna(k+2) == 3) %TAG
                tag = tag + 1;
            end
            
            if (dna(k) == 4 && dna(k+1) == 3 && dna(k+2) == 1) %TGA
                tga = tga + 1;
            end
            
            if (dna(k) == 4 && dna(k+1) == 1 && dna(k+2) == 1) %TAA
                taa = taa + 1;
            end
        
        end    
    end
end

% add the counting values of stop codons to the the array stop_freq
stop_freq = [stop_freq, tag];
stop_freq = [stop_freq, tga];
stop_freq = [stop_freq, taa];

% calculating average, min and max number of bases that being used
min_len = min (arr);
max_len = max (arr);
mean_len = mean (arr);


pro_seg = length (arr); % number of protein segments
percent = (sum(arr)/n)*100; % percentage of DNA being used in the protein process

%determine the min and max counting values of stop codons
most_stop = max (stop_freq);
least_stop = min (stop_freq);

%% Display the result
fprintf ("**DNA analysis**\n");
fprintf ("- Total protein-coding segments: %d \n", pro_seg);
fprintf ("- Average length: %5.2f \n", mean_len);
fprintf ("- Maximum length: %d \n", max_len);
fprintf ("- Minimum length: %d \n", min_len);
fprintf ("- %5.2f percent of this DNA is directly used in the protein coding process. \n\n",percent);


%% Determine which type of stop codon is the most used and the least used
switch (most_stop)
    case tag
        fprintf ("The most frequently used stop codon is TAG with %d uses.\n", most_stop);
    case tga
        fprintf ("The most frequently used stop codon is TGA with %d uses.\n", most_stop);
    case taa
        fprintf ("The most frequently used stop codon is TAA with %d uses.\n", most_stop);
    otherwise
end

switch (least_stop)
    case tag
        fprintf ("The least frequently used stop codon is TAG with %d uses.\n", least_stop);
    case tga
        fprintf ("The least frequently used stop codon is TGA with %d uses.\n", least_stop);
    case taa
        fprintf ("The least frequently used stop codon is TAA with %d uses.\n", least_stop);
    otherwise
end
