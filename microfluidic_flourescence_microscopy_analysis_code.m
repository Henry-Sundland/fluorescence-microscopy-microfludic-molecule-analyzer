%this code loops through tif files of molecule in channel
%tracks molecule via center of mass
%calculate mean squared displacement, diffusion, gyration tensor
%radius of gyration, molecular weight (bp) from calibration intensity....
%stores bunch of stuff into cell array, which can be seen at the bottom of
%the code
%code allows for user input, to make things faster for me...will ask you
%questions as it runs

%this is the first version of the code that I will post...I will clean it
%up and make it more user friendly soon, and will make a python
%version.....however, it isn't too difficult to use.
%any questions or suggestions, then hit me up! I love image analysis lol

%need to include the linefunc file and the linefitter file in the same
%folder as this code

intensity_calibration = input('Enter intensity calibration value to find molecular weight (should be in units of bp/intensity): ')
pixel_resolution = input('Enter pixel to nm conversion, based off of microscope resolution (pixel size in nm): ')



%reads values from csv files

myfiles_water = dir(''); %paste directory path to tif files here




filenames_water = {myfiles_water(:).name}';
tif_files_open = filenames_water(endsWith(filenames_water, '.tif'));
files_open = fullfile(tif_files_open);
%avg

h = 1;
movie_cell_array_raw_vids = {};
movie_cell_array_masked_vids = {};

% Creating strings in order to append to save to excel files
msd_y_string = 'msd_y_';
R_g_string = 'R_g_';
time_string = 'time_';
csv_for_string = '.csv';
net_velocity_for_string = 'net_v';
t_d_string = 't_d_';
t_c_string = 't_c_';
P_e_string = 'P_e_';
r_ratio = 'r_ratio_';
intensity_name = 'intensity_';
diffusion_string = 'diffusion_';
diffusion_uncertainty_string = 'diffusion_uncertainty_';
identity_file_string = 'identify_of_';
txt_string = '.txt';

masking_value = 5;



bool_vids = input('Enter 1 if you want option to watch vids: ')
bool_ex = 0;


while h <= size(tif_files_open,1)

               if bool_vids == 1
            
               bool_ex = input('Enter 1 if you want to play masked video: ' )
               disp('Now playing the current video in the data set: ')
                   
               end
            
               
               x = tif_files_open(h);
               y = char(x);
               filename = convertCharsToStrings(y);
               InfoImage=imfinfo(filename); 
                xdim=InfoImage(1).Width;  %xdim is the width of the image (in pixels)
                ydim=InfoImage(1).Height; %ydim is the height of the image (in pixels)
                frames=length(InfoImage); %how many frames
                movie=zeros(ydim,xdim,frames,'uint16');  %this is your 3D matrix of pixel intensities
                TifLink = Tiff(filename, 'r');
                colormap('gray')
                for i=1:frames
                    TifLink.setDirectory(i);
                    movie(:,:,i)=TifLink.read();
                end
                TifLink.close();
            
            
                movie_cell_array_raw_vids{h} = movie;
                frames_for_movie(h) = frames;
                filename
                for k = 1:frames
                    
                    %avg background test for frame
                    testing_movie = movie(:,:,k);
                    average_intense = mean(testing_movie);
                    max_intense = max(average_intense);
                    testing_movie(testing_movie >= max_intense) = 0;

                    % looping through frame and saving all intensities on
                    % each pixel that is not zero
                    t = 1;
                    intensities_no_zeroes = [];
                    dummy = 0;
                    for ugh = 1: size(testing_movie,1)
                        for ugh_2 = 1:size(testing_movie,2)
                            
                            dummy = testing_movie(ugh, ugh_2);

                            if dummy ~= 0
                                intensities_no_zeroes(t) = dummy;
                                dummy = 0;
                                t = t + 1;
                            end
                            dummy = 0;
                        end
                    end
                    dummy = 0;

                  background_avg = mean(intensities_no_zeroes);
%                     background_avg = mean(background_avg);
                    


              
                    % Noise subtraction of current movie iteration
                    img = movie(:,:,k) - background_avg;
                    movie_2 = img;
                    movie_3 = img;
                    
                    movie_3 = imgaussfilt(movie_2, 1);
                    average_intense = mean(movie_3);
                    max_intense = max(average_intense);
                    
                    movie_3(movie_3<max_intense + masking_value ) =0;
                    movie_3(movie_3 > 0) = 1;
                    
                    img = movie_2.*movie_3;
                    %masked_for_save(img(1),img(2),k);
            
                    %will play vid based on user input
                    if bool_ex == 1
                        imagesc(img)
                        pause(0.001)
            
                    end


          



            
                    
            
                    % Center of mass calculation
                    xcom = 0;
                    ycom = 0;
                    Itot=sum(img(:));
                    Intensity_size(k) = Itot; %for the size!
                    avg_intensity(k) = Itot; % to be used to find molecular weight of molecule
                    xInt = sum(img,1);
                    yInt = sum(img,2);
                    xcom = double(sum(xInt.*(1:(xdim)))./Itot);
                    ycom = double(sum(yInt.*(1:(ydim))')./Itot);
                    xloc(k) = xcom;
                    yloc(k) = ycom;
            
            
                    % gyration tensor, eigen values and R_g calculations
                    gxx = 0;
                    gyy = 0;
                    gxy = 0;
                    dummy = 0;
            
                    for i = 1:ydim
                        dummy = 0;
                        for j = 1:xdim
                            I = img(i,j);
                            dummy = double(I*(j-xcom)^2);
                            gxx = gxx + dummy;
                            dummy = double(I*(i-ycom)^2);
                            gyy = gyy + dummy;
                            dummy = double(I*(j-xcom)*(i - ycom));
                            gxy = gxy + dummy;
                            dummy = 0;
                        end
                    end
            
                    G(1,1) = gxx;
                    G(2,2) = gyy;
                    G(1,2) = gxy;
                    G(2,1) = gxy;
                    G = double(G/Itot); % Gyration tensor for first frame
                    %G_tot{k} = G;
                    
                    [V D] = eig(G); %eigen values D
                
                    first_eig = D(1,1);
                    second_eig = D(2,2);
                    R_g(k) = sqrt(first_eig + second_eig)*(172.4)*(1/1000);
            
            
            
                end


                 %prompting user if they wish to adjust intensity
                    bool_threshold = 0;

                    bool_threshold = input('Do you wish to adjust intensity threshold? Enter one for yes.')

                    if bool_threshold == 1
                        'current value is: ' 
                        
                        masking_value

                        masking_value = input('Enter new threshold integer value here:')

                    end

                
                %movie_cell_array_masked_vids{h} = masked_for_save;
            
                % converting pixels to micro meters
                for i = 1: frames
                    xloc(i) = xloc(i)*(pixel_resolution)*(1/1000);
                    yloc(i) = yloc(i)*(pixel_resolution)*(1/1000);
                end
                
                
                % converting frames to seconds
                
                for i = 1:frames
                    time_vec(i) = double(i*(1/50));
                end
            
            
                % % This here calculates and plots msd for x-dimension
                msd_x =zeros(1,frames);
                for j = 1:length(xloc)
                        for k = 1: (frames - j)
                            squarey_x = ((xloc(k+j) - xloc(k))^2)/(frames - j);
                            msd_x(j) = msd_x(j) + squarey_x;
                            
                        end
                        %store_x(j,:) = msd_x;
                end
            
            
            
               % This here calculates msd for y-dimension
                msd_y =zeros(1,frames);
                for j = 1:length(yloc)
                       for k = 1: (frames - j)
                           squarey_y = ((yloc(:,k+j) - yloc(:,k))^2)/(frames - j);
                           msd_y(j) = msd_y(j) + squarey_y;
                           
                       end
                       %store_y(j,:) = msd_y;
                end
            
            
            
                % finding velocities in x and y
                velocity_x = [];
                velocity_y = [];
                velocity_net = [];
                time_for_velocity = [];
                i = 1;
                
                while i <= frames - 1
                    velocity_x(i) = (xloc(i + 1) - xloc(i))/(time_vec(i + 1) - time_vec(i));
                    velocity_y(i) = (yloc(i + 1) - yloc(i))/(time_vec(i + 1) - time_vec(i));
                    velocity_net(i) = sqrt(velocity_y(i)^2 + velocity_x(i)^2);
                    i = i + 1;
                end
            
            
                average_R_g = mean(R_g);
            
    if bool_threshold ~= 1


                % Now finding diffusion for each msd(i); will do by fitting
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
            
                bool_plot = 0;
            
            
                bool_plot = input('Do you want to find the diffusion value? Enter 1 for yes: ')
            
             
            
            
            
                while bool_plot == 1
            
                    
                     'blue is msd_y and orange is msd_x'
                    'enter one to choose msd_y and enter two to choose msd_x; enter 3 to skip this one'
                    plot(msd_y)
                    hold on
                    plot(msd_x)
                    hold off
                    msd_choice = input('choose:')
                    
                    if msd_choice == 1
                        msd = msd_y;
                    end
                    
                    if msd_choice == 2
                        msd = msd_x;
                    end
            
          
            
            
            
            
                    plot(msd)



            
                    lower_limit = input('Enter lower limit integer to fit to: ')
                    upper_limit = input('Enter upper limit integer to fit to: ')
            
                    diffusion = msd(lower_limit:upper_limit)/(2*time_vec(lower_limit:upper_limit));
            
                   
                    P_e = 1;
                    % 
            
                    x_fit = time_vec(lower_limit:upper_limit);
                    y_fit = msd(lower_limit:upper_limit);
                    
           
                    
                    
                    w=ones(size(y_fit));
                    %w=1./dy.^2;
                    
                    %%this is where the meat happens
                    [fitout,resid,J,cov,mse]=nlinfit( x_fit,  y_fit, @linefunc, [1 1],'Weight',w);
                    ci = nlparci(fitout,resid,'jacobian',J);
                    
                    FitSlope=(ci(2,2)+ci(2,1))/2;
                    DeltaFitSlope=(ci(2,2)-ci(2,1))/2;
                    
                    FitInt=(ci(1,2)+ci(1,1))/2;
                    DeltaFitInt=(ci(2,2)-ci(1,1))/2;
                    
                    %% outputting and plotting
                    'slope'
                    [ FitSlope/2 DeltaFitSlope/2]
                    'intercept'
                    [ FitInt DeltaFitInt]
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    D = FitSlope/2;
                    D_uncertainty = DeltaFitSlope/2;
            
            
                   
                    
            
            
            
            
            
            
            
                    %don't worry about this stuff :)
                    diffusion_calc = mean(diffusion);
                    r_h = (k_boltz*T)/(6*pi*viscosity*D);
                    r_h = r_h*(1E+6);
                   
            
                    bool_plot = 0;
                end
            
            
            
            
            
            
            
            
            
            
            
            
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
              


                R = [];
                R = R_g/r_h;
            
            
               
            
                % converting intensity to bp
                
                base_pair = avg_intensity*(intensity_calibration);




                % writing to cell array
                C_data{1,h} = D; %diffusion for this row
                C_data{2,h} = D_uncertainty; %diffusion uncertainty
                C_data{3, h} = avg_intensity; % intensities of frames for this row
                C_data{4,h} = base_pair; %base pairs from intensities using lambda calibration
                C_data{5,h} = msd_x; %msdx for vid
                C_data{6,h} = msd_y; %msdy for vid
                C_data{7,h} = msd; % msd that was used to find diffusion for vid
                C_data{8,h} = R_g; %radius of gyration for vid
                C_data{9,h} = r_h; % hydrodynamic radius for vid....don't worry about this...yet.....
                C_data{10,h} = R; %ratio of r_g to r_h
                C_data{11,h} = time_vec; % time in seconds for vid
                C_data{12,h} = filename;% names for which vid that was used
             
            
                h = h + 1;

    end

end

save('molecule_data.mat', 'C_data');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%