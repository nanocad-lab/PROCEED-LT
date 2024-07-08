function found = check_hspice_conclude(filePath, targetLineFromEnd, targetString)
    % filePath: The file pathe of the file to be checked
    % targetLineFromEnd: The number of lines from the end of the file to search
    % targetString: The string to search for
    % found: found = 1 if the target string is found in the target line, otherwise found = 0

    fid = fopen(filePath, 'r');
    if fid == -1
        error('File can not be opened. Check the file path.');
    end
    
    % Locate the end of the file
    fseek(fid, 0, 'eof');
    fileSize = ftell(fid);
    
    linesFound = 0;
    position = fileSize;
    newlineChar = 10; % ASCII code for newline character

    % Search for the target line
    while linesFound <= targetLineFromEnd && position > 0
        position = position - 1;
        fseek(fid, position, 'bof');
        if fread(fid, 1, 'char') == newlineChar
            linesFound = linesFound + 1;
        end
    end

    % Go back to the target line position and read the target line
    fseek(fid, position + 1, 'bof');
    targetLine = fgets(fid);

    fclose(fid);
    
    % Check if the target string is found in the target line
    found = contains(targetLine, targetString);
end
