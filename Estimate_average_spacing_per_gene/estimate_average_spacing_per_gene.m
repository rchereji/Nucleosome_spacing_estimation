function estimate_average_spacing_per_gene(dyads_filename)
% This function takes the MATLAB filename that contains the Dyads
% distribution (e.g. created using by the get_dyads_and_occupancy.m
% function) and estimates the average nucleosome spacing near the 5' end of
% each gene
%
% Inputs:
% dyads_filename - MATLAB file containing the dyads distribution in a cell 
%                  array called Dyads
%
% Output: This function saves a file called
% sprintf('Classification_%s_%d_bp.mat', data_set, a), which contains the
% following info: 
% ORF      - List of the ORFs
% Spacing  - estimated nucleosome spacing corresponding to each ORF
% Shift    - the estimated shift of the nucleosome array
% bestCorr - a score that gives the correlation between the nucleosome
%            array and the test oscillatory profile that best fitted the data
%
% Example:
% estimate_average_spacing_per_gene('Dyads_WT_A_120_160.mat') - estimate
% the spacing in WT cells, replicate 1

beforeRef = 200;
afterRef = 1000;
AnnotationsFile = 'Yeast_Annotations_NAR_2016.mat';

% Load the dyads distribution
load(dyads_filename, 'Dyads');

% Name of the data set
data_set = dyads_filename(7 : end-4);

% Extend all dyads to a footprint of 101 bp
a = 101;
Occ = extend_dyads_to_footprint(Dyads, a);

% Normalize the coverage by the chromosome average
noChr = numel(Occ);
for chr = 1:noChr
    tmp = Occ{chr};
    
    % Mark the regions with very high coverage (sequencing artifacts)
    Filter = tmp > 100 * nanmean(tmp);
    % Mark the rDNA region (region with very high coverage, which skews the chromosome average)
    if chr == 12
        Filter(451000 : 469000) = true;
    end
    % Dilate the filter
    Filter = imdilate(Filter, strel('line', 101, 0));
    
    % Rescale the coverage by the average of the non-filtered/"good" regions
    tmp(Filter) = nan;
    Occ{chr} = Occ{chr} / nanmean(tmp);
end

% Chromosome lengths
chrLen = [230218,813184,316620,1531933,576874,270161,1090940,562643,...
    439888,745751,666816,1078177,924431,784333,1091291,948066];
% chrName = {'chrI';'chrII';'chrIII';'chrIV';'chrV';'chrVI';'chrVII';'chrVIII';...
%     'chrIX';'chrX';'chrXI';'chrXII';'chrXIII';'chrXIV';'chrXV';'chrXVI'};

% Align +1 nucleosomes for all genes
load(AnnotationsFile, 'ORF', 'Chr', 'Plus1', 'TTS', 'Watson');

% Keep only the genes with available annotations for the typical +1 nuc. pos.
idx = (~isnan(Plus1));
ORF = ORF(idx);
Chr = Chr(idx);
Plus1 = Plus1(idx);
TTS = TTS(idx);
Watson = Watson(idx);

noGenes = numel(ORF);
AlignedProfile = nan(noGenes, 1 + beforeRef + afterRef);
for g = 1:noGenes
    if Watson(g)
        leftEdge = max([Plus1(g) - beforeRef, 1]);
        rightEdge = min([Plus1(g) + afterRef, chrLen(Chr(g))]);
        AlignedProfile(g, beforeRef + 1 - (Plus1(g) - leftEdge)...
            : beforeRef + 1 + (rightEdge - Plus1(g))) = ...
            Occ{1,Chr(g)}(leftEdge : rightEdge);
    else
        leftEdge = max([Plus1(g) - afterRef, 1]);
        rightEdge = min([Plus1(g) + beforeRef, chrLen(Chr(g))]);
        AlignedProfile(g, beforeRef + 1 - (rightEdge - Plus1(g))...
            : beforeRef + 1 + (Plus1(g) - leftEdge)) = ...
            fliplr(Occ{1,Chr(g)}(leftEdge : rightEdge));
    end
end

goodGenes = mean(AlignedProfile, 2) > 0.4;
AlignedProfile = AlignedProfile(goodGenes, :);
ORF = ORF(goodGenes);
Chr = Chr(goodGenes);
Plus1 = Plus1(goodGenes);
TTS = TTS(goodGenes);
Watson = Watson(goodGenes);
DistToTTS = (TTS - Plus1) .* (2 * Watson - 1);

% Eliminate genes for which the TTS is too close to +1 nucleosome and less
% than 2 nucleosomes can fit before TTS
goodGenes = DistToTTS > 200;
AlignedProfile = AlignedProfile(goodGenes, :);
ORF = ORF(goodGenes);
DistToTTS = DistToTTS(goodGenes);

% Use at most 800 bp (5 nucs.) or 600 bp (4 nucs.) to estimate the spacing
DistToTTS = min(DistToTTS, 600);

% Generate oscillatory pattern made of 10 Gaussians, with sigma=40 and
% separated by a distance d (variable)
Pattern = cell(1, 220);
for d = 130:220
    Pattern{d} = gaussmf(-beforeRef:afterRef,[40 0]) + ...
        gaussmf(-beforeRef:afterRef,[40 d]) + ...
        gaussmf(-beforeRef:afterRef,[40 2*d]) + ...
        gaussmf(-beforeRef:afterRef,[40 3*d]) + ...
        gaussmf(-beforeRef:afterRef,[40 4*d]) + ...
        gaussmf(-beforeRef:afterRef,[40 5*d]) + ...
        gaussmf(-beforeRef:afterRef,[40 6*d]) + ...
        gaussmf(-beforeRef:afterRef,[40 7*d]) + ...
        gaussmf(-beforeRef:afterRef,[40 8*d]) + ...
        gaussmf(-beforeRef:afterRef,[40 9*d]);
end

% Estimate the spacing for each gene
noGenes = numel(ORF);
Shift = nan(noGenes, 1);
Spacing = nan(noGenes, 1);
bestCorr = nan(noGenes, 1);

parfor g = 1:noGenes
    bestCorrTest = -1;
    Profile = AlignedProfile(g, :);
    bpToRemove = afterRef - DistToTTS;
    if bpToRemove > 0
        Profile(end - bpToRemove : end) = NaN; % Discard the region downstream of TTS
    end
    
    for d = 130:220 % spacing to test
        for shift = -70:0 % Pattern shifted to the left
            R = corrcoef([Profile(1:end+shift)', Pattern{d}(1-shift:end)'],'rows','pairwise'); 
            C =  R(1,2);
            if C > bestCorrTest
                bestCorrTest = C;
                Shift(g) = shift;
                Spacing(g) = d;
            end
        end
        
        for shift = 1:60 % Pattern shifted to the right
            R = corrcoef([Profile(1+shift:end)', Pattern{d}(1:end-shift)'],'rows','pairwise'); 
            C =  R(1,2);
            if C > bestCorrTest
                bestCorrTest = C;
                Shift(g) = shift;
                Spacing(g) = d;
            end
        end
    end
    
    bestCorr(g) = bestCorrTest;
%     disp(g) % Uncomment to display the genes that were already processed
end

save(sprintf('Classification_%s_%d_bp.mat', data_set, a), 'ORF', 'Spacing', 'Shift', 'bestCorr')
