<script>
  let bg = -1;

  function getBg() {
    fetch("./diabetes")
      .then(d => d.json())
    //   .then(d => (rand = d));
      .then(function (d) {
		  console.log(d)
		  return bg = d
	  });
  }


// Copied from internet!
  	import { onMount } from "svelte";
	import Chart from "chart.js";
	onMount(async () => {});
	function renderChart(bg) {
		var ctx = document.getElementById("myChart").getContext("2d");
		let x = [...Array(10).keys()]
		var chart = new Chart(ctx, {
		type: "line",
		data: {
			// labels: [0, 1, 2, 3, 4, 5, 6],
			labels: x, 
			datasets: [
			{
				label: "My First dataset",
				backgroundColor: "rgb(255, 99, 132)",
				borderColor: "rgb(255, 99, 132)",
				fill: 'False',
				data: bg
				// data: x
			}
			]
		},
		options: {}
		});

  }

// End copied from internet!
</script>

<h1>Your number is {}!</h1>
<button on:click={getBg}>Get a random number</button>

<button on:click={renderChart(bg)}>Load</button> 
<canvas id="myChart"></canvas> 
