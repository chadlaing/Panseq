function add() {
				  transfer($("#source-list li.ui-selected"));
				}
				function addAll() {
				  transfer($("#source-list li:visible"));
				}
				function transfer(listItems) {
				  listItems.fadeOut(function() {
				    $(this)
				      .removeClass("ui-selected")
				      .clone()
				      .appendTo("#target-list")
				      .fadeIn()
				      .data("index", $("#source-list li").index($(this)))
				      .highlight();
				  });
				}



function remove() {
  $("#target-list li.ui-selected").fadeOut(function() {
    $("#source-list li")
      .eq($(this).data("index"))
      .removeClass("ui-selected")
      .fadeIn()
      .highlight();
 
     $(this).remove();
  });
}


function removeAll() {
  $("#target-list li").fadeOut()
    .promise().done(function() {
      $("#target-list li").remove();
      $("#source-list li:hidden").fadeIn().highlight();
    });
}


function addHiglightPlugin() {
  $.fn.highlight = function() {
    return this
      .addClass("li-transfer-highlight")
      .removeClass("li-transfer-highlight", 400);
  }
}

