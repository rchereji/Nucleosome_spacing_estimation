function [frag_length, frag_count] = get_DNA_fragment_lengths(bam_filename)
% This function computes the read length distribution from the BAM file
%
% Inputs:
% bam_filename - BAM file containing the genomic alignments of the
% sequenced DNA fragments
%
% Output:
% frag_length - length of a DNA fragment, between 1 and 1000 bp
% frag_count  - the number of reads corresponding to each DNA size
%
% Example:
% [frag_length, frag_count] = get_DNA_fragment_lengths('WT_A.bam');
%

% Check the inputs
if nargin == 0
    error('You didn''t provide the data file!')
else
    if (exist(bam_filename, 'file') ~= 2)
        error('File "%s" does not exist in the current folder!', bam_filename)
    end
end

chrLen = [230218,813184,316620,1531933,576874,270161,1090940,562643,...
    439888,745751,666816,1078177,924431,784333,1091291,948066];
chrName = {'chrI';'chrII';'chrIII';'chrIV';'chrV';'chrVI';'chrVII';'chrVIII';...
    'chrIX';'chrX';'chrXI';'chrXII';'chrXIII';'chrXIV';'chrXV';'chrXVI'};
noChr = numel(chrName);
all_fragment_lengths = [];

% Process all chromosomes
for chr = 1 : noChr
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
    
    all_fragment_lengths = [all_fragment_lengths; fragmentLengths];
    
    % Stop when at least 100,000 proper paired-end reads have been analyzed
    % If all the fragments need to be considered, just comment out the
    % following IF statement
    if numel(all_fragment_lengths > 100000)
        break
    end
end

frag_length = [1:1000]';
frag_count = arrayfun(@(z) sum(ismember(all_fragment_lengths, z)), frag_length);

setName = bam_filename(1:end-4);
save(sprintf('Length_histogram_%s.mat', setName), 'frag_length', 'frag_count');

% Plot the figure
figure('Position', [50, 50, 400, 300]);
plot(frag_length, 100*frag_count/sum(frag_count), 'LineWidth', 1);
xlim([0 500])
set(gca, 'XTick', 0:100:500, 'XGrid', 'on', 'GridLineStyle', '--')
set(gca, 'XMinorTick', 'on', 'TickLength', [0.03 0.025])
xlabel('Fragment length (bp)')
ylabel('Percentage (%)')
title({'Length histogram'; sprintf('Sample: %s', setName)}, 'interpreter', 'none')
box off

% Create inset with mono-nuc data
axes('Position',[.5 .43 .35 .32])
plot(frag_length, 100*frag_count/sum(frag_count), 'LineWidth', 2);
xlim([100 200])
set(gca, 'XTick', 100:20:200, 'XGrid', 'on', 'GridLineStyle', '--')
title({'Zoom-in view:'; '100 bp - 200 bp'})

% Save the figure as an EPS file
set(gcf, 'PaperPositionMode', 'auto')
print(gcf, sprintf('Length_histogram_%s', setName), '-depsc', '-painters');
