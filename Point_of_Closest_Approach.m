function PL_intersection = Point_of_Closest_Approach(WL_list, PL_1, PL_2)
    x1 = WL_list;
    y1 = PL_1;
    x2 = WL_list;
    y2 = PL_2;
    
    dist_min = 100;
    idx_min_dist = [1, 1];
    for i = 1:length(x1)
        for j = 1:length(x2)
            euclid_dist = sqrt((y2(j) - y1(i))^2 + (x2(j) - x1(i))^2);
            if euclid_dist < dist_min
                dist_min = euclid_dist;
                idx_min_dist = [i, j];
            end
        end
    end
    
    PL_intersection = (PL_1(idx_min_dist(1)) + PL_2(idx_min_dist(2))) / 2;
    
end
