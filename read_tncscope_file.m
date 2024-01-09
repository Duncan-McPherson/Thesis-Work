function [data] = read_tncscope_file(file)
%READ_TNCSCOPE_FILE Reads TNCScope data (*.sco)
%   Detailed explanation goes here

line = '';
NCh = 0; 
 
%split up file part
[filepath,name,ext] = fileparts(file);

%assert, that this is really a scope file
assert(strcmp(ext, '.sco'));

% open file
[fID,errmsg] = fopen(file);

%only take action, if fopen was successful
if fID > -1
    
    %% read input file
    
    % read all additional information
    while(~strcmp(line, '%[DATA]'))
        line = fgetl(fID);
        
        % global settings
        if contains(line,'%[GLOBAL]')
            while(~is_line_empty(line)) %global description ends at empty line
                line = fgetl(fID); 
                
                %globalinfo =
            end
        end
        
        % channels
        if contains(line,'%[CURVE')
            %new channel found --> count
            NCh = NCh +1;
            
            %loop over all channel lines
            while(~is_line_empty(line)) %channel description ends at empty line
                line = fgetl(fID);
                
                % save the information on the channel
                if startsWith(line, '%CHANNELTYPE')            
                    channelinfo(NCh).type = extractAfter(line, '=');
                elseif startsWith(line, '%AXISTYPE')
                    channelinfo(NCh).axis = extractAfter(line, '=');
                elseif startsWith(line, '%DIMENSION')
                    channelinfo(NCh).unit = extractAfter(line, '=');
                elseif startsWith(line, '%COUNT')
                    channelinfo(NCh).Ns = str2double(extractAfter(line, '='));
                elseif startsWith(line, '%CYCLETIME')
                    channelinfo(NCh).Ts = str2double(extractAfter(line, '='))*1e-6; % in  [s]
                end
            end
        end
    end
    
    % set up format string for textscan depending on the number of found
    % channels
    formatstring = '';
    for i=1:NCh+1 % channels + time data
        formatstring = [formatstring, '%f;'];
    end
    
    % read all time and measurement data
    filedata = textscan(fID, formatstring);
    
    % close file
    fclose(fID);    
    
    
    %% construct output data     
    data.filename = strcat(name,ext);
    data.path = filepath;
    
    for i=1:NCh
        data.channel(i).time = [0:1:channelinfo(i).Ns-1] * channelinfo(i).Ts;
        data.channel(i).values = filedata{i+1}(1:channelinfo(i).Ns);
        
        %output row vectors only
        data.channel(i).time = data.channel(i).time(:);
        data.channel(i).values = data.channel(i).values(:);
        
        
        data.channel(i).name = get_channel_name(channelinfo(i));
        data.channel(i).unit = get_channel_unit(channelinfo(i));
        data.channel(i).axis = channelinfo(i).axis;
        
        data.channel(i).Ts = channelinfo(i).Ts;
        data.channel(i).fs = 1/channelinfo(i).Ts;
        data.channel(i).Ns = channelinfo(i).Ns;
    end
    
    
else
    disp(errmsg);
end


%% help functions
    function [bReturn] = is_line_empty(line)
        bReturn = strcmp(strtrim(line), '%');
    end
    
    function [bReturn] = is_debug_channel(channeltype)
        bReturn = strcmp(strtrim(channeltype), 'DSP-Debug');
    end

    function [name] = get_channel_name(channelinfo)
        if is_debug_channel(channelinfo.type)
            switch (channelinfo.unit)
                case 'id(10)'
                    name = 'Torque Generating Current';
                case 'id(222)'
                    name = 'Commanded position';
                case 'id(217)'
                    name = 'Actual position linear encoder';
                case 'id(125)'
                    name = 'Commanded velocity';
                case 'id(126)'
                    name = 'Actual velocity motor encoder';
            end
        else
            name = channelinfo.type;
        end
    end

    function [unit] = get_channel_unit(channelinfo)
        if is_debug_channel(channelinfo.type)
            switch (channelinfo.unit)
                case 'id(10)'
                    unit = 'A';
                case 'id(222)'
                    unit = 'int';
                case 'id(217)'
                    unit = 'int';
                case 'id(125)'
                    unit = 'rad/sec';
                case 'id(126)'
                    unit = 'rad/sec';

            end
        else
            unit = channelinfo.unit;
        end
    end

%% print data information
print_TNC_channel_info(data);

end

