graph_ts_200ka <- function(display = TRUE, datasets = c("data_odp_967_22", "data_odp721_722_terr", "data_odp_709", "data_icdp_chb", "data_kl09", "data_kl11", "data_kl15", "data_lake_tana"), index = 1) {
    d <- get(datasets[index])

    #*******************************200ka Plot of wentness index*************************************
    # Create the plot
    png(filename = paste0("graphs/200ka_polygonplot", datasets[index],".png"))

    yt <- d[,2] # your time series y values
    xt <- d[,1] # your time series age values

    # Assuming d is a matrix with two columns
    mean_yt <- mean(d[!is.na(d[,2]), 2], na.rm = TRUE)
    std_dev <- sd(d[!is.na(d[,2]), 2], na.rm = TRUE)

    # Define ylim based on 3-sigma range
    y_min <- mean_yt - 5 * std_dev
    y_max <- mean_yt + 5 * std_dev
    
    #***************************Plotting the graph********************************
    if (!display) {
        plot(xt, yt, type = "n", xlab = "", ylab = "", xlim = range(xt), ylim = range(yt), main = paste0("200ka Plot of ", datasets[index]))
    }

    # Create the filled areas
    # age matrix
    fill_x <- c(xt, rev(xt))
    # Terrigenous matrix
    fill_y <- c(yt, rep(mean_yt, length(yt)))
    fill_y[yt < mean_yt] <- mean_yt
    fill_y[1] <- mean_yt
    fill_y[length(fill_y)] <- mean_yt
    color_rgb <- col2rgb('#D95319') / 255
    polygon(fill_x, fill_y, col = rgb(color_rgb[1], color_rgb[2], color_rgb[3], alpha = 0.1), border = NA)

    fill_y2 <- c(yt, rep(mean_yt, length(yt)))
    fill_y2[yt > mean_yt] <- mean_yt
    fill_y2[1] <- mean_yt
    fill_y2[length(fill_y2)] <- mean_yt
    color_rgb <- col2rgb('#0072BD') / 255
    polygon(fill_x, fill_y2, col = rgb(color_rgb[1], color_rgb[2], color_rgb[3], alpha = 0.1), border = NA)

        # Close the graphics device
    dev.off()
    
    # #*******************alternative ggplot graph********************************
    # library(ggplot2)

    # # Assuming d is a data frame with two columns
    # mean_val <- mean(d[!is.na(d[,2]), 2], na.rm = TRUE)
    # std_dev <- sd(d[!is.na(d[,2]), 2], na.rm = TRUE)

    # # Define ylim based on 3-sigma range
    # y_min <- mean_val - 5 * std_dev
    # y_max <- mean_val + 5 * std_dev

    # # Create a data frame for the filled areas
    # fill_pos_data <- data.frame(x = c(d[,1], rev(d[,1])), 
    #                         y = c(d[,2], rep(mean_val, nrow(d))))
    # # negative values are replaced with the mean value
    # fill_pos_data$y[fill_data$y < mean_val] <- mean_val

    # # Create a data frame for the filled areas
    # fill_neg_data <- data.frame(x = c(d[,1], rev(d[,1])), 
    #                         y = c(d[,2], rep(mean_val, nrow(d))))
    # # positive values are replaced with the mean value
    # fill_neg_data$y[fill_data$y >= mean_val] <- mean_val


    # # Create the plot
    # p <- ggplot(d, aes(x = yt, y = xt)) +
    # geom_line() +
    # geom_polygon(data = fill_pos_data, aes(x = x, y = y), 
    #             fill = rgb(217/255, 83/255, 25/255, alpha = 0.9), 
    #             color = NA) +
    # geom_polygon(data = fill_neg_data, aes(x = x, y = y), 
    #             fill = rgb(173/255, 216/255, 230/255, alpha = 0.9), 
    #             color = NA) +
    # ylim(y_min, y_max) +
    # theme_minimal() +
    # labs(x = "Time (kyrs BP)", y = "", title = paste0("200ka Plot of Wentness Index of ", datasets[1]))

    # ggsave(filename = paste0("200ka_plot_",datasets[index],".png"), plot = p, path = "graphs/", width = 10, height = 6, dpi = 300)
}
