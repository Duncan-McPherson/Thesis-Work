function [] = print_TNC_channel_info(data)
%PRINT_TNC_CHANNEL_INFO Prints read channel information of TNCOpt and
%TNCScope files

fprintf('Read TNC file %s\n', fullfile(data.path,data.filename))
for i=1:length(data.channel)
    fprintf('Channel %i:\n',i);
    fprintf('\tName: %s\n', data.channel(i).name);
    fprintf('\tAxis: %s\n', data.channel(i).axis);
    fprintf('\tUnit: %s\n', data.channel(i).unit);
    fprintf('\t  Ts: %.2e\n', data.channel(i).Ts);
    fprintf('\t  Ns: %i\n\n', data.channel(i).Ns);
end

end

