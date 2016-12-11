function [Data,TimeStamp,DataSmooth] = generateData(M,N)
    d = 5;  % size of the object
    X = 1:N; 
    Y = 1:M;
    AER_TH = 20;    
    TimeStamp =[];
    
    for i=1:1:N-d
        
        x = X-i;y=Y;
        x = mod(x,N);
        y = mod(y,M);
        %countourX = (heaviside(x-d)-heaviside(x)).*(x).*(x-d-1);
        %countourY = (heaviside(y-d)-heaviside(y)).*(y).*(y-d-1);
        countourX = 10*(heaviside(x-d)-heaviside(x)).*ones(size(x));
        countourY = 10*(heaviside(y-d)-heaviside(y+1)).*ones(size(y));

        frame  = countourY'*countourX;
        if(i>1)
            frame = frame+wgn(M,N,1);
            gradient = frame-frame_;
            gradientD = sign(gradient.*floor(abs(gradient/AER_TH)));
            RawImgInput(:,:,i) = frame;
            DataSmooth(:,:,i) = gradient;
            Data(:,:,i)= gradientD;
            TimeStamp = [TimeStamp i];
            image(gradient,'CDataMapping','scaled');
        end

        frame_ = frame;
        pause(0.01)
    end
    
    
    