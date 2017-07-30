var m = 'PA LUT model normalization parameters:' +
    'Amp_pa %%  1.1428861359Offset_pa 0.0000000000 '
if(m.indexOf("%%")!=-1){
    console.log(m.replace(/%%/g,'='));
}else{
    console.log("none");
}


function rea(a,b){
    console.log(a);
    if(b){
        console.log(b);
    }
    return;
}
rea(3);
rea(3,9);
rea(4,6,10);