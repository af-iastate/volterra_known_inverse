function [sfdr,sfdrBounds] = calcSFDR(x,spacing,frqs,bounds,sRate,plt)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% fix this later
if isempty(sRate)
    fs = 125e6;
else
    fs = sRate;
end

% get spacing numbers
offset = spacing(1);
duration = spacing(2);
skip = spacing(3);

% assumes cell array
% nFiles = length(files);

% number of (two-tone) signals
nPulses = size(frqs,1);

if nargin<6
    plt = zeros(nFiles,1);
elseif length(plt)==1
    plt = repmat(plt,nPulses,1);
end

% label bins with frequencies
fd = fs/duration;

% find frequency bins for bounds
lbf = round(bounds(1)/fd)+1;
ubf = round(bounds(2)/fd)+1;

% set up matrices
%nmse = zeros(nFiles,nPulses);
sfdr = zeros(nFiles,nPulses);
sfdrBounds = zeros(nFiles,nPulses);

nh = ceil(duration/2);

% for n = 1:nFiles
%     % load new file
%     x = input_adc_samples_4ch(1,files{n});
    
    % loop over pulses
    for h = 1:nPulses
        % pull out segment
        iStart = (h-1)*skip+offset;
        iEnd = iStart + duration - 1;
        xc = db(fft(x(iStart:iEnd)));
        
        % remove dc and fs/2
        %xc(1:lbf) = -inf; % covered below
        %xc(ubf+1:nh-1) = -inf;
        xc(1) = -inf;
        xc = xc(1:nh-1);
        
        % do stats
        % determine bin of tones
        cfr = frqs(h,:);
        cfr(cfr==0) = []; % eliminate dc; encodes no tone
        cfr = round(cfr./fd)+1;
        
        [minTonePow, minToneIdx] = min(xc(cfr));
        minToneIdx = cfr(minToneIdx);
        xc(cfr) = -inf;
        [maxSpurPow, maxSpurIdx] = max(xc);
        [maxSpurBounds, maxBSpurIdx] = max(xc(lbf+1:ubf));
        maxBSpurIdx = maxBSpurIdx + lbf;
        
        % fill in matrix
        sfdr(n,h) = minTonePow - maxSpurPow;
        sfdrBounds(n,h) = minTonePow - maxSpurBounds;
        
        % optional plotting
        if plt(h)
            xc = db(fft(x(iStart:iEnd)));
            fr = 0:fd:(length(xc)-1)*fd;
            figure; plot(fr,xc)
            hold all
            plot(fr(minToneIdx),xc(minToneIdx),'*')
            plot(fr([maxSpurIdx,maxBSpurIdx]),xc([maxSpurIdx,maxBSpurIdx]),'o')
            ax = axis;
            plot(fr(repmat(lbf+1,3,1)),[ax(3),(ax(3)+ax(4))/2,ax(4)],'k-')
            plot(fr(repmat(ubf,3,1)),[ax(3),(ax(3)+ax(4))/2,ax(4)],'k-')
            
            % process title
            tl = files{1};
            if length(tl)>40
                tl = [tl(1:17),' ... ',tl(end-17:end)];
            end
            title(sprintf('%s, Burst: %d',tl,h))
        end
        
    end
    
% end


end

