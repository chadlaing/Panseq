$(document).ready(function() {
   
   // the hoverable tabs with jQueryUI
   $(function() {
    $( "#tabs" ).tabs({
      event: "click"
    });
  });

   // submit button with jQueryUI
   $(function() {
    $( "input[type=submit], a, button" )
      .button()
      .click(function( event ) {
        event.preventDefault();
      });
  });

  //actually submit when submitButton clicked
  $("#formSubmit").click(function() {
        submitForm();
  });

  //follow nav links
  $("#hrefHome").button().click(function() {
      window.location.href='/panseq/home/';
  });

  $("#hrefAnalyses").button().click(function() {
      window.location.href='/panseq/analyses/';
  });

  $("#hrefContact").button().click(function() {
      window.location.href='/panseq/contact/';
  });

  //Panseq github link
  $("#hrefPanseqGithub").button().click(function() {
      window.open('https://github.com/chadlaing/Panseq');
  });

  //Panseq BMC Bioinformatics link
  $("#hrefPanseqBMC").button().click(function() {
    window.open('http://www.biomedcentral.com/1471-2105/11/461');
  });
  

   // init for the list transfer
  $("#query-source-list, #query-target-list").selectable();
  $("#query-add-button").click(queryAdd);
  $("#query-remove-button").click(queryRemove);
  $("#query-remove-all-button").click(queryRemoveAll);

  $("#reference-source-list, #reference-target-list").selectable();
  $("#reference-add-button").click(referenceAdd);
  $("#reference-remove-button").click(referenceRemove);
  $("#reference-remove-all-button").click(referenceRemoveAll);

  //selected styling events for list items
 $("#query-source-list, #query-target-list, #reference-source-list, #reference-target-list").selectable({
  unselected: function(){
    $(":not(.ui-selected)", this).each(function(){
      $(this).removeClass('ui-state-highlight');
    });
  },
  selected: function(){
    $(".ui-selected", this).each(function(){
      $(this).addClass('ui-state-highlight');
    });
  }
});

  //hover styling
  $(".listContainer ul li, #query-target-list, #reference-target-list").hover(
      function () {
        $(this).addClass('ui-state-hover');
      }, 
      function () {
        $(this).removeClass('ui-state-hover');
      }
  );


  //populate hidden list on form submit

}); //end document ready

  //submit novel form
  function submitForm(){

  //get all li items from the query-target-list
  //iterate over all of them
  //store the name attribute from the list, not the displayed text
  //join the name attributes together and add to the hidden field
  //this hidden field is what is passed from the form
   var queryArray=[];
   $("#query-target-list li").each(function(){  
       queryArray.push($(this).attr('name'));
    });  
    var queryList = queryArray.join();
    $('input[name=querySelected]').val(queryList);

    //reference strains
    var refArray=[];
   $("#reference-target-list li").each(function(){  
       refArray.push($(this).attr('name'));
    });  
    var refList = refArray.join();
    $('input[name=referenceSelected]').val(refList);

    this.submit();
  }


  // //add and add-all buttons
  function queryAdd() {
    queryTransfer($("#query-source-list li.ui-selected"));
  }

  function queryTransfer(listItems) {
    listItems.fadeOut(0,function() {
      $(this)
        .removeClass("ui-state-highlight ui-selected")
        .clone()
        .appendTo("#query-target-list")
        .fadeIn(0)
        .data("index", $("#query-source-list li").index($(this)));
    });
  }

  //remove button
  function queryRemove() {
  // $("#target-list li.ui-selected").fadeOut(function() {
    $("#query-target-list li.ui-selected").fadeOut(0,function() {
    $("#query-source-list li")
      .eq($(this).data("index"))
      .removeClass("ui-selected")
      .fadeIn(0);    
     $(this).remove();
    });
  }


  //remove all
  function queryRemoveAll() {
  $("#query-target-list li").fadeOut(0)
    .promise().done(function() {
      $("#query-target-list li").remove();
      $("#query-source-list li:hidden").fadeIn(0);
    });
  }

  // //add and add-all buttons
  function referenceAdd() {
    referenceTransfer($("#reference-source-list li.ui-selected"));
  }

  function referenceTransfer(listItems) {
    listItems.fadeOut(0,function() {
      $(this)
        .removeClass("ui-state-highlight ui-selected")
        .clone()
        .appendTo("#reference-target-list")
        .fadeIn(0)
        .data("index", $("#reference-source-list li").index($(this)));
    });
  }

  //remove button
  function referenceRemove() {
  // $("#target-list li.ui-selected").fadeOut(function() {
    $("#reference-target-list li.ui-selected").fadeOut(0,function() {
    $("#reference-source-list li")
      .eq($(this).data("index"))
      .removeClass("ui-selected")
      .fadeIn(0);    
     $(this).remove();
    });
  }

  //remove all
  function referenceRemoveAll() {
  $("#reference-target-list li").fadeOut(0)
    .promise().done(function() {
      $("#reference-target-list li").remove();
      $("#reference-source-list li:hidden").fadeIn(0);
    });
  }


