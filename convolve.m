[input, fs] = wavread('voice.wav');
ir = wavread('ir.wav');

inputlen = length(input);
irlen = length(ir);
outputlen = inputlen+irlen-1;

N = 1024;
ol = 4;
olchunk = N/ol;
hw = hann(N);
numblocks = ceil(inputlen/olchunk);

% add zeros to the end of the input signal until we reach a multiple of the block size
inputpadded = zeros(numblocks*N, 1);
inputpadded(1:inputlen) = input;

irspec = fft(ir);
outputbuffer = zeros(outputlen, 1);

for i=1:numblocks,

	start = (i-1)*olchunk+1;
	finish = start+N-1;
	
	inputbuf = zeros(irlen, 1);
	inputbuf(1:N) = inputpadded(start:finish) .* hw;
	
	inputspec = fft(inputbuf);
	
	specprod = inputspec .* irspec;
	
	output = ifft(specprod);
	
	outputbuffer(start:start+irlen-1) = outputbuffer(start:start+irlen-1) .+ output;

end

% normalize
divisor = max(abs(outputbuffer));
divisor = 1.0/divisor; % use recip to avoid divides below
outputbuffer = outputbuffer*divisor;