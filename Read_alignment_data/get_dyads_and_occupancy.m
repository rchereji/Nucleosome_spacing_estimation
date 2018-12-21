function [Dyads, Occ] = get_dyads_and_occupancy(bam_filename, Lmin, Lmax)
% This function loads all the alignments from a BAM file, and computes the
% coverage profile (Occ) and the distribution of the fragment centers (Dyads)
%
% Inputs:
% bam_filename - BAM file containing the genomic alignments of the
% sequenced DNA fragments
% Lmin / Lmax  - DNA sizes to be considered (e.g. 120-160 bp)
% 
% Output: 
% Dyads - cell array (one element for each chromosome) containing the
%         number of nucleosome centers (dyads) that were detected at each
%         genomic position 
% Occ   - cell array containing the number of DNA fragments that cover
%         each genomic position (i.e. occupancy / coverage profile)
%
% Example:
% [Dyads, Occ] = get_dyads_and_occupancy('WT_A.bam', 120, 160);
%

% Check the inputs
if nargin == 0
    error('You didn''t provide the data file!')
end

if nargin < 3
    if (exist(bam_filename, 'file') ~= 2)
        error('File "%s" does not exist in the current folder!', bam_filename)
    else
        % Default size range (no size selection)
        Lmin = 120;
        Lmax = 160;
    end
end

% Initialize the Occ & Dyads cell arrays
genomeVer = 'sacCer3';
chrLen = [230218,813184,316620,1531933,576874,270161,1090940,562643,...
    439888,745751,666816,1078177,924431,784333,1091291,948066];
chrName = {'chrI';'chrII';'chrIII';'chrIV';'chrV';'chrVI';'chrVII';'chrVIII';...
    'chrIX';'chrX';'chrXI';'chrXII';'chrXIII';'chrXIV';'chrXV';'chrXVI'};
noChr = numel(chrName);
Occ = cell(1, noChr);
Dyads = cell(1, noChr);

fprintf('Starting to import data from "%s".\n', bam_filename)
fprintf('Size selection: %d <= L <= %d\n', Lmin, Lmax)

% Process all chromosomes
NoReadsPerChr = zeros(1, noChr);
parfor chr = 1 : noChr
    % Create BioMap object
    bm = BioMap(bam_filename, 'SelectReference', chrName{chr});
    
    % Eliminate the reads with very low quality (QMAP = 0)
    bm = getSubset(bm, getMappingQuality(bm) > 0);
    
    % Filter out the reads which fall outside the reference
    bm = getSubset(bm, getStop(bm) <= chrLen(chr));
    
    % Get the reads that map to the Watson strand
    Indices = filterByFlag(bm, 'pairedInMap', true, 'strandQuery', 0);
    bm_filtered_Watson = getSubset(bm, Indices);
    
    % Get the reads that map to the Crick strand
    Indices = filterByFlag(bm, 'pairedInMap', true, 'strandQuery', 1);
    bm_filtered_Crick = getSubset(bm, Indices);
    
    % Match the pairs (1 Watson read + 1 Crick read) that form a paired-end read
	[~, Watson_Idx, Crick_Idx] = intersect(bm_filtered_Watson.Header, bm_filtered_Crick.Header);
    
    % Compute the fragment lengths
    leftBP = getStart(bm_filtered_Watson, Watson_Idx);
    rightBP = getStop(bm_filtered_Crick, Crick_Idx);
    fragmentLengths = rightBP - leftBP + 1;
    
    % Size selection, acording to the specified limits (amin, amax)
    goodInd = ((fragmentLengths >= Lmin) & (fragmentLengths <= Lmax));
    leftBP = leftBP(goodInd);
    rightBP = rightBP(goodInd);
    
    noGoodPairs = numel(leftBP);
    
    % Construct Dyads distribution
    Centers = round((rightBP + leftBP)/2);
    uniqueCenters = unique(Centers);
    [~, Index] = ismember(Centers, uniqueCenters);
    NumberUniqueCenter = histc(Index, 1:numel(uniqueCenters));
    
    Dyads{chr} = zeros(1, chrLen(chr));
    Dyads{chr}(uniqueCenters) = NumberUniqueCenter;
    
    % Compute Occ
    uniqueLeftBP = unique(leftBP);
    [~, Index] = ismember(leftBP, uniqueLeftBP);
    NumberUniqueleftBP = histc(Index, 1:numel(uniqueLeftBP));
    
    OccDerivative = zeros(1, chrLen(chr) + 1); % add 1 position for the case when some reads have the right end exactly at the end of chr
    OccDerivative(uniqueLeftBP) = NumberUniqueleftBP;
    
    uniqueRightBP = unique(rightBP);
    [~, Index] = ismember(rightBP, uniqueRightBP);
    NumberUniqueRightBP = histc(Index, 1:numel(uniqueRightBP));
    
    OccDerivative(uniqueRightBP + 1) = OccDerivative(uniqueRightBP + 1) - NumberUniqueRightBP';
    tmp = cumsum(OccDerivative);
    
    Occ{chr} = tmp(1 : end-1);
    
    NoReadsPerChr(chr) = noGoodPairs;
    fprintf('Chr. %d done.\n', chr);
end

% Create a log string
infoStr = sprintf('File: %s\nReference    Number of Pairs    Density (pairs/bp)\n', bam_filename);

for chr = 1 : noChr
    % Fill in the log string
    infoStr = sprintf('%s%8s\t\t%d\t\t\t\t%0.6f \n', infoStr,...
        chrName{chr}, NoReadsPerChr(chr), double(NoReadsPerChr(chr) / chrLen(chr)));
end 
TotalNoReads = sum(NoReadsPerChr);
infoStr = sprintf('%s%8s\t\t%d\t\t\t\t%0.6f \n', infoStr,...
    'Total', TotalNoReads, double(TotalNoReads / sum(chrLen)));

% Save the Occupancy and Dyads cell arrays
setName = bam_filename(1:end-4);
save(sprintf('Occupancy_%s_%d_%d.mat', setName, Lmin, Lmax), 'Occ', 'TotalNoReads', 'genomeVer');
save(sprintf('Dyads_%s_%d_%d.mat', setName, Lmin, Lmax), 'Dyads', 'TotalNoReads', 'genomeVer');

% Save the log file with the sequencing depth
fileID = fopen(sprintf('%s_log_%d_%d.txt', bam_filename, Lmin, Lmax),'w');
fprintf(fileID,'%s', infoStr);
fclose(fileID);

fprintf('File "%s" has been successfully processed.\n\n', bam_filename)
