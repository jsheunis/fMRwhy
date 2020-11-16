
param_bids_dataset;

function setupReportStructure() {

  allrunsKeys = [];
  allrunsVals = [];
  allrunsSessions = []; 
  allrunsTasks = [];
  allrunsRuns = [];

  bids_dataset['report_runs'].forEach(function(sestaskrun, index) {

    temp_key = sestaskrun;

    var session, task, run;
    var s = sestaskrun.indexOf('ses-');
    if (s >= 0) {
      session = sestaskrun.slice(s+4).split('_')[0];
    } else {
      session = '';
    };
    var t = sestaskrun.indexOf('task-')
    if (t >= 0) {
      task = sestaskrun.slice(t+5).split('_')[0]
    } else {
      task = '';
    };
    var r = sestaskrun.indexOf('run-')
    if (r >= 0) {
      run = sestaskrun.slice(r+4).split('_')[0]
    } else {
      run = '';
    };

    allrunsKeys.push(temp_key);
    allrunsVals.push(temp_key);
    allrunsSessions.push(session);
    allrunsTasks.push(task);
    allrunsRuns.push(run);
  });
  // console.log(allrunsKeys)
  // console.log(allrunsVals)
  // If map_rois==false, remove nodes/elements related to ROIs
  if (!bids_dataset['map_rois']) {
    var element = document.getElementById('anatroitext');
    element.remove();
    var element2 = document.getElementById('anatroiimgs');
    // var node= document.getElementById("parent");
    element2.querySelectorAll('*').forEach(n => n.remove());
    element2.remove();
  }

  // If map_rois==false, remove nodes/elements related to ROIs
  if (!bids_dataset['include_physio']) {
    document.getElementById('physqcplots').remove();
    document.getElementById('physqcplots_img').remove();
  }

  // Create options for runselector from allruns array, and create objects to access sessions/tasks/runs from
  allrunsObj = {};
  allrunsSessionsObj = {};
  allrunsTasksObj = {};
  allrunsRunsObj = {};
  var sel = document.getElementById('runselector');
  allrunsVals.forEach(function(txt, index) {
    // create new option element
    var opt = document.createElement('option');
    // create text node to add to option element (opt)
    opt.appendChild(document.createTextNode(txt) );
    // set value property of opt
    opt.value = allrunsKeys[index]; 
    // add opt to end of select box (sel)
    sel.appendChild(opt);
    allrunsObj[allrunsKeys[index]] = txt;
    allrunsSessionsObj[allrunsKeys[index]] = allrunsSessions[index];
    allrunsTasksObj[allrunsKeys[index]] = allrunsTasks[index];
    allrunsRunsObj[allrunsKeys[index]] = allrunsRuns[index];
  });

  // Stats table setup
  var newRow, newCell, newText;
  var tableRef = document.getElementById('statsTable');

  allrunsKeys.forEach(function(key, index) {
    newRow = tableRef.insertRow(-1);
    newCell = newRow.insertCell(-1);
    newText = document.createTextNode(key);
    newCell.appendChild(newText);

    console.log(bids_dataset['all_runs_stats'][index])

    Object.keys(bids_dataset['all_runs_stats'][index]).forEach(function(statKey, statIndex) {
      newCell = newRow.insertCell(-1);
      newText = document.createTextNode(bids_dataset['all_runs_stats'][index][statKey]);
      newCell.appendChild(newText);
      newCell.align = "right";
    });
  });
  
}


$(document).ready(function() {
  $('#runselector').change(function(){


    var imgArrayExt = {
      "mean_img": "_mean.png",
      "std_img": "_std.png",
      "tsnr_img": "_tsnr.png",
      // "grayplot_ro_img": "_echo-2_desc-RO_grayplot.png",
      // "grayplot_gso_img": "_echo-2_desc-GSO_grayplot.png",
      // "grayplot_lmotor_img": "_echo-2_desc-leftMotor_grayplot.png",
      // "grayplot_biamy_img": "_echo-2_desc-bilateralAmygdala_grayplot.png"
    };

    var grayplotsDesc = {
      "grayplot_ro_img": "RO",
      "grayplot_gso_img": "GSO",
      // "grayplot_lmotor_img": "_echo-2_desc-leftMotor_grayplot.png",
      // "grayplot_biamy_img": "_echo-2_desc-bilateralAmygdala_grayplot.png"
    };

    if (bids_dataset['map_rois']) {
      bids_dataset['roi_desc'].forEach(function(roi_desc, index) {
        grayplots["grayplot_" + roi_desc + "_img"] = roi_desc;
      });
    }

    // Reset sources for images
    // $("#" + key).attr("src", "img/sub-" + bids_dataset['sub'] + '_' + $(this).val() + imgArray[key]);
    for (var key in imgArrayExt) {

      fn = fmrwhy_bidsjs_constructFilename({
        filetype:'func',
        sub:bids_dataset['sub'],
        ses:allrunsSessionsObj[$(this).val()],
        task:allrunsTasksObj[$(this).val()],
        acq:'', rec:'',
        run:allrunsRunsObj[$(this).val()],
        echo:'',
        space:'individual',
        desc:'',
        ext:imgArrayExt[key]});
      $("#" + key).attr("src", "img/" + fn);
    }

    var echo = '';
    if (bids_dataset['is_multiecho']) {
      echo = bids_dataset['template_echo']
    }
    for (var key in grayplotsDesc) {
      fn = fmrwhy_bidsjs_constructFilename({
        filetype:'func',
        sub:bids_dataset['sub'],
        ses:allrunsSessionsObj[$(this).val()],
        task:allrunsTasksObj[$(this).val()],
        acq:'', rec:'',
        run:allrunsRunsObj[$(this).val()],
        echo:echo,
        space:'', 
        desc:grayplotsDesc[key],
        ext:'_grayplot.png'});

      $("#" + key).attr("src", "img/" + fn);
    }



    $("#physqcplots_img").attr("src", "img/" + bids_dataset['sub'] + '_' + $(this).val() + bids_dataset['physio_str']);
    // Reset heading names
    $("#spatialqcplots").html( "Spatial QC Plots: " + allrunsObj[$(this).val()] );
    $("#tempqcplots").html("Temporal QC plots: " + allrunsObj[$(this).val()] );
    $("#physqcplots").html("PhysIO QC plots: " + allrunsObj[$(this).val()] );

  });
});


function fmrwhy_bidsjs_constructFilename({filetype="", sub="", ses="", task="", acq="", rec="", run="", echo="", space="", desc="", ext=""}={}) {
  // BIDS 1.4.0
  // This does not yet take into account whether an entity is OPTIONAL or REQUIRED for a specific file modality/type (e.g. 'func', 'anat', 'meg', etc)
  // It assumes all parameters are optional and just constructs the filename given the input
  // also, file extension is passed as parameter ==> should be automatically determined from file modality/type
  //{'ses', 'acq', 'rec', 'run', 'echo'} ==> to keep in line with fMRwhy implementation

  var filetypes = ['anat', 'func']; // currently supported in fMRwhy
  var entities = ['sub', 'ses', 'task', 'acq', 'rec', 'run', 'echo', 'space', 'desc', 'ext'];

  var filename = '';
  var newStr;
  entities.forEach(function(key, index) {
    if (!(eval(key)==="")) {
      if ((key=='sub') || (key=='ext')) {
        newStr = eval(key);
      } else {
        newStr = '_' + key + '-' + eval(key);
      }
      filename += newStr;
    }
  });
  
  return filename;
}


function openTab1(evt, tabName) {
  // Declare all variables
  var i, tabcontent, tablinks;

  // Get all elements with class="tabcontent" and hide them
  tabcontent = document.getElementsByClassName("tabcontent1");
  for (i = 0; i < tabcontent.length; i++) {
    tabcontent[i].style.display = "none";
  }

  // Get all elements with class="tablinks" and remove the class "active"
  tablinks = document.getElementsByClassName("tablinks1");
  for (i = 0; i < tablinks.length; i++) {
    tablinks[i].className = tablinks[i].className.replace(" active", "");
  }

  // Show the current tab, and add an "active" class to the button that opened the tab
  document.getElementById(tabName).style.display = "block";
  evt.currentTarget.className += " active";
}


function openTab2(evt, tabName) {
  // Declare all variables
  var i, tabcontent, tablinks;

  // Get all elements with class="tabcontent" and hide them
  tabcontent = document.getElementsByClassName("tabcontent2");
  for (i = 0; i < tabcontent.length; i++) {
    tabcontent[i].style.display = "none";
  }

  // Get all elements with class="tablinks" and remove the class "active"
  tablinks = document.getElementsByClassName("tablinks2");
  for (i = 0; i < tablinks.length; i++) {
    tablinks[i].className = tablinks[i].className.replace(" active", "");
  }

  // Show the current tab, and add an "active" class to the button that opened the tab
  document.getElementById(tabName).style.display = "block";
  evt.currentTarget.className += " active";
}

function openTab3(evt, tabName) {
  // Declare all variables
  var i, tabcontent, tablinks;

  // Get all elements with class="tabcontent" and hide them
  tabcontent = document.getElementsByClassName("tabcontent3");
  for (i = 0; i < tabcontent.length; i++) {
    tabcontent[i].style.display = "none";
  }

  // Get all elements with class="tablinks" and remove the class "active"
  tablinks = document.getElementsByClassName("tablinks3");
  for (i = 0; i < tablinks.length; i++) {
    tablinks[i].className = tablinks[i].className.replace(" active", "");
  }

  // Show the current tab, and add an "active" class to the button that opened the tab
  document.getElementById(tabName).style.display = "block";
  evt.currentTarget.className += " active";
}

function openTab4(evt, tabName) {
  // Declare all variables
  var i, tabcontent, tablinks;

  // Get all elements with class="tabcontent" and hide them
  tabcontent = document.getElementsByClassName("tabcontent4");
  for (i = 0; i < tabcontent.length; i++) {
    tabcontent[i].style.display = "none";
  }

  // Get all elements with class="tablinks" and remove the class "active"
  tablinks = document.getElementsByClassName("tablinks4");
  for (i = 0; i < tablinks.length; i++) {
    tablinks[i].className = tablinks[i].className.replace(" active", "");
  }

  // Show the current tab, and add an "active" class to the button that opened the tab
  document.getElementById(tabName).style.display = "block";
  evt.currentTarget.className += " active";
}