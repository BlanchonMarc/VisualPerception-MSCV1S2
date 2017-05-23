%loadtestimages    load test images
%    [M,nomia,nomib]=loadtestimages(image)
%
% by X. Armangue
% (c) Mr3D - University of Girona, September 2002
%
function [M,nomia,nomib]=loadtestimages(image)

switch image,
    case 1,
        fitxer='imatges/calibnoni.dat';
        nomia='imatges/calibnonia.tif';
        nomib='imatges/calibnonib.tif';
    case 2,
        fitxer='imatges/inria.dat';
        nomia='imatges/inriaa.tif';
        nomib='imatges/inriab.tif';
        colormap('gray');
    case 3,
        fitxer='imatges/j.dat';
        nomia='imatges/j1.tif';
        nomib='imatges/j2.tif';
    case 4,
        fitxer='imatges/jordif.dat';
        nomia='imatges/jordifa.tif';
        nomib='imatges/jordifb.tif';
    case 5,
        fitxer='imatges/palamos.dat';
        nomia='imatges/palamosa.tif';
        nomib='imatges/palamosb.tif';
    case 6,
        fitxer='imatges/paret.dat';
        nomia='imatges/pareta.tif';
        nomib='imatges/paretb.tif';
    case 7,
        fitxer='imatges/patro3d.dat';
        nomia='imatges/patro3da.bmp';
        nomib='imatges/patro3db.bmp';
    case 8,
        fitxer='imatges/piri.dat';
        nomia='imatges/piria.tif';
        nomib='imatges/pirib.tif';
    case 9,
        fitxer='imatges/robot.dat';
        nomia='imatges/robota.tif';
        nomib='imatges/robotb.tif';
    case 10,
        fitxer='punts/punts00o00.dat';
        nomia=[];
        nomib=[];
    case 11,
        fitxer='punts/punts00o05.dat';
        nomia=[];
        nomib=[];
    case 12,
        fitxer='punts/punts00o10.dat';
        nomia=[];
        nomib=[];
    case 13,
        fitxer='punts/punts01o00.dat';
        nomia=[];
        nomib=[];
    case 14,
        fitxer='punts/punts01o05.dat';
        nomia=[];
        nomib=[];
    case 15,
        fitxer='punts/punts01o10.dat';
        nomia=[];
        nomib=[];
    case 16,
        fitxer='punts/punts05o00.dat';
        nomia=[];
        nomib=[];
    case 17,
        fitxer='punts/punts05o05.dat';
        nomia=[];
        nomib=[];
    case 18,
        fitxer='punts/punts05o10.dat';
        nomia=[];
        nomib=[];
    case 19,
        fitxer='punts/punts10o00.dat';
        nomia=[];
        nomib=[];
    case 20,
        fitxer='punts/punts10o05.dat';
        nomia=[];
        nomib=[];
    case 21,
        fitxer='punts/punts10o10.dat';
        nomia=[];
        nomib=[];
    case 22,
        fitxer='punts/newpunts00o00.dat';
        nomia=[];
        nomib=[];
    case 23,
        fitxer='punts/newpunts00o05.dat';
        nomia=[];
        nomib=[];
    case 24,
        fitxer='punts/newpunts00o10.dat';
        nomia=[];
        nomib=[];
    case 25,
        fitxer='punts/newpunts01o00.dat';
        nomia=[];
        nomib=[];
    case 26,
        fitxer='punts/newpunts01o05.dat';
        nomia=[];
        nomib=[];
    case 27,
        fitxer='punts/newpunts01o10.dat';
        nomia=[];
        nomib=[];
    case 28,
        fitxer='punts/newpunts05o00.dat';
        nomia=[];
        nomib=[];
    case 29,
        fitxer='punts/newpunts05o05.dat';
        nomia=[];
        nomib=[];
    case 30,
        fitxer='punts/newpunts05o10.dat';
        nomia=[];
        nomib=[];
    case 31,
        fitxer='punts/newpunts10o00.dat';
        nomia=[];
        nomib=[];
    case 32,
        fitxer='punts/newpunts10o05.dat';
        nomia=[];
        nomib=[];
    case 33,
        fitxer='punts/newpunts10o10.dat';
        nomia=[];
        nomib=[];
end
fid = fopen(fitxer,'r');
M = fscanf(fid,'%f',[4,inf]);
fclose(fid);
disp(sprintf('Points loaded from file %s',fitxer));
