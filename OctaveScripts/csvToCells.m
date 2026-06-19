function [c] = csvToCells(fileName)

  c = {};

  if ~exist(fileName, 'file')
    fprintf('File %s not found.', fileName);
    c = -1;
    return
  end

  f1 = fopen(fileName, 'r');

  eofFound = false;

  currentLine = 0;

  while ~eofFound

    strIn = fgetl(f1);

    if numel(strIn) == 1
       if (strIn(1, 1) - 0) == -1

         eofFound = true;

      end
    end

    if ~eofFound

      currentLine = currentLine + 1;

      idCommas = find(strIn == ',');

       id1 = [1, (idCommas + 1)];
       id2 = [(idCommas - 1), size(strIn, 2)];

       for currentColumn = 1:size(id1, 2)

         c{currentLine, currentColumn} = strIn(1, id1(1, currentColumn):id2(1, currentColumn));

       end

    end

  end

  fclose(f1);

end
