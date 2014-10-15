#!/usr/bin/env Rscript

merge.rows <- function(r1, r2) {
    if (is.na(r1["start"])) return(list(next.seg=r2, new.rows=NULL))

    # make sure rows are in correct order
    if ((r2[1] < r1[1]) | (r2[1] == r1[1] & r2[2] < r1[2])) {
        cat("swapping rows\n")
        tmp <- r2
        r2 <- r1
        r1 <- tmp
    }

    d1 <- ifelse(is.na(r1[3]), r2[3], r1[3])
    d2 <- ifelse(is.na(r1[4]), r2[4], r1[4])
    
    # intervals don't start at the same point
    if (r1[1] < r2[1]) {
        # case 1: intervals don't overlap
        if (r1[2] < r2[1]) {
            new.rows <- r1
            next.seg <- r2
        }
    
        # case 2: r1 ends where r2 starts
        # case 3: r1 ends in the middle of r2
        else if (r1[2] < r2[2]) {
            new.rows <- rbind(c(r1[1], r2[1]-1, r1[3:4]),
                              c(r2[1], r1[2], d1, d2))
            next.seg <- c(r1[2]+1, r2[2], r2[3:4])
        }
    
        # case 4: r1 ends at the same point as r2
        else if (r1[2] == r2[2]) {
            new.rows <- rbind(c(r1[1], r2[1]-1, r1[3:4]),
                              c(r2[1], r2[2], d1, d2))
            next.seg <- rep(NA, 4)
        }
    
        # case 5: r1 ends after r2
        else if (r1[2] > r2[2]) {
            new.rows <- rbind(c(r1[1], r2[1]-1, r1[3:4]),
                              c(r2[1:2], d1, d2))
            next.seg <- c(r2[2]+1, r1[2], r1[3:4])
        } else {
            cat("something is wrong\n", file=stderr())
            quit(status=1)
        }
    } else if (r1[1] == r2[1]) {
        # case 6: same start, r1 ends before r2
        if (r1[2] < r2[2]) {
            new.rows <- data.frame(r1[1:2], d1, d2)
            next.seg <- c(r1[2]+1, r2[2:4])
        }

        # case 7: same start, same end
        else if (r1[2] == r2[2]) {
            new.rows <- data.frame(r1[1:2], d1, d2)
            next.seg <- rep(NA, 4)
        } else {
            cat("something is wrong\n", file=stderr())
            quit(status=1)
        }
    }

    colnames(new.rows) <- c("start", "end", colnames(r1)[3:4])
    list(new.rows=new.rows, next.seg=next.seg)
}

merge.segs <- function (d1, d2) {
    d2[,colnames(d1)[3]] <- NA
    d1[,colnames(d2)[3]] <- NA
    d <- rbind(d1, d2)
    d <- d[order(d$start, d$end),]
    new.d <- do.call(rbind, lapply(1:(nrow(d)-1), function (i) {
        res <- merge.rows(d[i,], d[i+1,])
        d[i+1,] <<- res[["next.seg"]]
        res[["new.rows"]]
    }))
    if (is.na(tail(d, 1)[1]))
        as.data.frame(new.d)
    else
        as.data.frame(rbind(new.d, tail(d, 1)))
}
